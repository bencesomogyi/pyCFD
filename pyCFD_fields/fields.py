"""
module for abstract variable fields
"""
__author__ = "Bence Somogyi"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Bence Somogyi"
__email__ = "bencesomogyi@ivt.tugraz.at"
__status__ = "Prototype"

import pyCFD_config.config as Config_
import numpy
import sys
import pyCFD_mesh.patch as patch
if sys.platform == 'win32':
    import pyCFD_operators.cython_boost_win32.cy_operators as cy_operators
#    import pyCFD_general.cython_boost_win32.cy_general as cy_general
elif sys.platform == 'linux2':
    import pyCFD_operators.cython_boost_linux2.cy_operators as cy_operators
#    import pyCFD_general.cython_boost_linux2.cy_general as cy_general
else:
    sys.exit("unknown platform " + sys.platform + " found in geomTools.py, stopping...")
 
# "global variables"
#import __builtin__
__FIELDDIR__ = Config_.__FIELDDIR__

class Field:
    """
    basic class for fields
    """
    def __init__(self):
        self.name = ""
        """string for field name"""
        self.type = ""
        """string for field type: scalar or vector"""
        self.father = []
        """reference to owner mesh object"""
        self.patches = []
        """references to field patches"""
        self.A = []
        """surface values"""
            
    def get_patch(self, patch_name):
        """
        return a reference to the field patch with matching name
        
        :param patch_name: patch name of geometric patch (without '__fieldName')
        :type patch_name: str
        """
        patch_found = -1
        for patch_i,patch_ in enumerate(self.patches):
            if (patch_name+"__"+self.name) == patch_.name:
                patch_found = patch_i
                break
        
        if patch_found == -1:
            print patch_name + " does not exist in the mesh, possible patch names are:"
            for patch_ in self.patches:
                print patch_.name
            sys.exit()
        
        return self.patches[patch_found]
            
    def dot_Sf(self):
        """
        return the dot products of face vector values with face normals as
        scalar vector
        """
        if self.type != "vector":
            sys.exit(".dot_Sf() applied on non vector field!")
        face_number = len(self.father[0].faces)
        dot_product = numpy.zeros((face_number, 1))
        for i_ in range(face_number):
            dot_product[i_] = numpy.dot(self.A[i_], self.father[0].faces[i_].Sf)
        return dot_product
            
    def dot_abs_Sf(self):
        """
        return the dot products of face vector values with absolute face
        normals as scalar vector
        """
        if self.type != "vector":
            sys.exit(".dot_Sf() applied on non vector field!")
        face_number = len(self.father[0].faces)
        dot_product = numpy.zeros((face_number, 1))
        for i_ in range(face_number):
            dot_product[i_] = numpy.dot(self.A[i_], numpy.absolute(self.father[0].faces[i_].Sf))
        return dot_product
                        
class VolumeField(Field):
    """
    class for volume fields
    """
    def __init__(self):
        Field.__init__(self)
        self.V = []
        """cell values"""
    
    def update_boundary_values(self):
        """
        update boundary values
        """
        scalar_field = True
        if self.type == "vector":
            scalar_field = False

        mesh_ = self.father[0]
        patch_face_ids = []
        for i_ in xrange(len(mesh_.patches)):
            patch_face_ids.append(mesh_.patches[i_].ids)
        for face_ in mesh_.faces:
            if face_.isBnd == False:
                continue
            patch_face_i = face_.inPatchId
            patch_name = mesh_.patches[face_.bndId].name
            field_patch = self.get_patch(patch_name)
            if field_patch.type == "fixedValue":
                self.A[face_.id] = field_patch.values[patch_face_i]
            else: # "fixedGradient"
                vect_cell_to_face = numpy.add(face_.C,-face_.cells[0].C)
                # dist_cell_to_face = numpy.linalg.norm(vect_cell_to_face)
                dist_cell_to_face = cy_operators.cy_linalg_norm(vect_cell_to_face)
                # face_unit_vector = face_.Sf / numpy.linalg.norm(face_.Sf)
                face_unit_vector = face_.Sf / cy_operators.cy_linalg_norm(face_.Sf)
                
                field_grad = self.get_patch(patch_name).values[patch_face_i]
                cell_id = face_.cells[0].id
                delta_ = numpy.dot(field_grad * face_unit_vector, vect_cell_to_face) * dist_cell_to_face
                if scalar_field:
                    self.A[face_.id] = self.V[cell_id] + delta_
                else:
                    for component_ in range(3):
                        self.A[face_.id][component_] = self.V[cell_id][component_] + delta_
                        
    def load_init_fields(self, name_=""):
        """
        load initial field values and boundary conditions from saved .npy files
        
        :param name_: name of directory to load the files from located in :const:`pyCFD_config.config.__FIELDDIR__`
        :type name_:  str
        """
        if name_ == "":
            self.V = numpy.load(__FIELDDIR__+self.name+".npy")
        else:
            self.V = numpy.load(__FIELDDIR__+name_+"/"+self.name+".npy")
        for patch_ in self.patches:
            if name_ == "":
                patch_.values = numpy.load(__FIELDDIR__+patch_.name+".npy")
            else:
                patch_.values = numpy.load(__FIELDDIR__+name_+"/"+patch_.name+".npy")
        self.update_boundary_values()
        
    def copy_field(self, other_field):
        self.V[:] = other_field.V[:]
        self.A[:] = other_field.A[:]
        
        
class ScalarField(VolumeField):
    """
    class for scalar fields
    """
    def __init__(self, mesh_, field_name, init_value=0.):
        """
        **Constructor**
        
        :param mesh_: owner mesh object
        :type mesh_:  :class:`pyCFD_mesh.generic_mesh.GenericMesh`
        :param field_name: name of resulting field
        :type field_name: str
        :param init_value: default: 0, initial value
        :type init_value: float
        """
        VolumeField.__init__(self)
        
        self.name = field_name
        self.type = "scalar"
        
        self.V = numpy.ones(len(mesh_.cells)) * init_value
        self.A = numpy.ones(len(mesh_.faces)) * init_value
        
        self.father.append(mesh_)
        
        self.patches = []
        """reference to patches of field"""
        
        # initialize empty patches
        for patch_ in mesh_.patches:
            temp_field_patch = patch.FieldPatch(self, patch_)
            self.patches.append(temp_field_patch)
            
    def initialize_cell_with_vector(self, values):
        """
        set field cell values to the values given as input
        """
        if len(self.V) != len(values):
            print "can not initialize field with different sized vector!"
            print "field length: " + str(len(self.V))
            print "init length: " + str(len(values))
            sys.exit()
            
        for cell_i,field_value in enumerate(values):
            self.V[cell_i] = field_value
            
    def initialize_face_with_vector(self, values):
        """
        set field face values to the values given as input
        """
        if len(self.A) != len(values):
            print "can not initialize field with different sized vector!"
            print "field length: " + str(len(self.A))
            print "init length: " + str(len(values))
            sys.exit()
            
        for face_i,field_value in enumerate(values):
            self.A[face_i] = field_value
            

class VectorField(VolumeField):
    """
    class for vector fields
    """
    def __init__(self, mesh_, field_name, init_value=numpy.zeros(3)):
        """
        **Constructor**
        
        :param mesh_: owner mesh object
        :type mesh_:  :class:`pyCFD_mesh.generic_mesh.GenericMesh`
        :param field_name: name of resulting field
        :type field_name: str
        :param init_value: default: numpy.zeros(3), initial value
        :type init_value: numpy.array
        """
        VolumeField.__init__(self)
        
        self.name = field_name
        self.type = 'vector'
        
        if len(init_value) != 3:
            sys.exit("vectorField init value must be 3 element vector!")
            
        self.V = numpy.ones((len(mesh_.cells), 3))
        self.A = numpy.ones((len(mesh_.faces), 3))
        for i in range(3):
            self.V[:, i] *= init_value[i]
            self.A[:, i] *= init_value[i]
            
        self.father.append(mesh_)
        
        self.patches = []
        """reference to patches of field"""
        
        # initialize empty patches
        for patch_ in mesh_.patches:
            temp_field_patch = patch.FieldPatch(self, patch_)
            self.patches.append(temp_field_patch)
            
    def get_component_as_scalar_field(self, component_):
        """
        return one component of a volume vector field as a volume scalar field
        """
        names = ["_X", "_Y", "_Z"]
        field_component = ScalarField(self.father[0], self.name+names[component_])
        field_component.initialize_cell_with_vector(self.V[:, component_])
        field_component.initialize_face_with_vector(self.A[:, component_])
        for patch_i,patch_ in enumerate(self.father[0].patches):
            if self.patches[patch_i].type == "fixedValue":
                field_component.get_patch(patch_.name).set_patch_distributed(self.patches[patch_i].values[:, component_], self.patches[patch_i].type)
            else:
                field_component.get_patch(patch_.name).set_patch_distributed(self.patches[patch_i].values, self.patches[patch_i].type)
            
        return field_component
                
class SurfaceScalarField(Field):
    """
    class for surface scalar fields
    """
    def __init__(self, mesh_, field_name, init_value=0.0):
        """
        **Constructor**
        
        :param mesh_: owner mesh object
        :type mesh_:  :class:`pyCFD_mesh.generic_mesh.GenericMesh`
        :param field_name: name of resulting field
        :type field_name: str
        :param init_value: default: 0, initial value
        :type init_value: float
        """
        Field.__init__(self)
        
        self.name = field_name
        self.type = 'scalar'
        self.A = numpy.ones(len(mesh_.faces)) * init_value
        self.father.append(mesh_)
        
    def initialize_with_vector(self, values):
        """
        set field face values to the values given as input
        """
        if len(self.A) != len(values):
            print "can not initialize field with different sized vector!"
            print "field length: " + str(len(self.V))
            print "init length: " + str(len(values))
            sys.exit()
            
        for face_i,field_value in enumerate(values):
            self.A[face_i] = field_value
            
    def copy_field(self, other_field):
        self.A[:] = other_field.A[:]
        
class SurfaceVectorField(Field):
    """
    class for surface vector fields
    """
    def __init__(self, mesh_, field_name, init_value=numpy.zeros(3)):
        """
        **Constructor**
        
        :param mesh_: owner mesh object
        :type mesh_:  :class:`pyCFD_mesh.generic_mesh.GenericMesh`
        :param field_name: name of resulting field
        :type field_name: str
        :param init_value: default: numpy.zeros(3), initial value
        :type init_value: numpy.array
        """      
        Field.__init__(self)
        
        self.name = field_name
        self.type = 'vector'
        self.A = []
        
        if len(init_value) != 3:
            sys.exit("vectorField init value must be 3 element vector!")
            
        self.A = numpy.ones((len(mesh_.faces), 3))
        for i in range(3):
            for j in range(len(mesh_.faces)):
                self.A[j, i] *= init_value[i]
            
        self.father.append(mesh_)
        
    def initialize_with_vector(self, values):
        """
        set field face values to the values given as input
        """
        if len(self.A) != len(values):
            print "can not initialize field with different sized vector!"
            print "field length: " + str(len(self.V))
            print "init length: " + str(len(values))
            sys.exit()
            
        for face_i,field_value in enumerate(values):
            for component in range(3):
                self.A[face_i][component] = field_value[component]
            
    def copy_field(self, other_field):
        self.A[:] = other_field.A[:]
