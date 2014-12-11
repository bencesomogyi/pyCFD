"""
module for patches
"""
__author__ = "Bence Somogyi"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Bence Somogyi"
__email__ = "bencesomogyi@ivt.tugraz.at"
__status__ = "Prototype"

import sys
import numpy

class Patch:
    """
    class for geometric patches
    """
    def __init__(self, face_list, patch_name, check_independent=False):
        """
        **Constructor**
        
        :param face_list:         list of faces in the patch
        :type face_list:          :class:`pyCFD_mesh.face.Face`
        :param patch_name:        name to use for the resulting patch
        :type patch_name:         string
        :param check_independent: default: False, whether to check if faces are independent
        :type check_independent:  bool
        """
        if check_independent:
            facesIndependent = True
            for indFaceI in range(len(face_list)-1):
                for faceJ in face_list[indFaceI+1:]:
                    if faceJ.are_faces_equal(face_list[indFaceI],faceJ):
                        facesIndependent = False
            if facesIndependent == False:
                print "faces are not independent in ", self
                return
        
        self.faces = face_list
        """list of faces in the patch"""
        self.name = patch_name
        """name of patch"""
        self.ids = [face_list[face_i].id for face_i in xrange(len(face_list))]
        """list of face ids in the patch"""
                    
class FieldPatch:
    """class for field patches"""
    def __init__(self, field, geometric_patch):
        """
        **Constructor**
        
        The field patch is initialized with fixedValue type and with face
        values of 0.0.        
        
        :param field:           reference to field object
        :type field:            :class:`pyCFD_fields.fields.VolumeField`
        :param geometric_patch: reference to geometric patch
        :type geometric_patch:  :class:`Patch`
        """
        self.father = []
        """reference to father geometric patch object"""
        self.type = "fixedValue"
        """patch type: 'fixedValue' or 'fixedGradient'"""
        if field.type == "scalar":
            self.values = numpy.zeros(len(geometric_patch.faces))
        else:
            self.values = numpy.zeros((len(geometric_patch.faces),3))
        """list of boundary values for each patch faces"""
        self.field = field
        """reference to owner field"""
        
        self.father.append(geometric_patch)
        self.name = geometric_patch.name + "__" + field.name
    
    def set_patch_uniform(self, boundary_value, boundary_type):
        """
        set up patch with uniform boundary values
        """
        # check type
        available_types = []
        available_types.append("fixedValue")
        available_types.append("fixedGradient")
        
        if boundary_type not in available_types:
            print "not supported boundary type "+boundary_type+" defined, availabe types are:"
            print available_types
            sys.exit()

        self.type = boundary_type
            
        vector_field = False
        if (self.field.type == "vector"):
            vector_field = True
            
        fixed_value = False
        if boundary_type == "fixedValue":
            fixed_value = True
        
        if vector_field and fixed_value:
            self.values = numpy.zeros((len(self.father[0].faces),3))
            for face_i,face_ in enumerate(self.father[0].faces):
                self.values[face_i, :] = boundary_value
                face_.bndType = boundary_type
        else:
            self.values = numpy.zeros(len(self.father[0].faces))
            for face_i,face_ in enumerate(self.father[0].faces):
                self.values[face_i] = boundary_value
                face_.bndType = boundary_type
                        
                    
    def set_patch_distributed(self, boundary_values, boundary_type):
        """
        set up patch with non-uniform boundary values
        """
        # check boundary type
        if boundary_type != "fixedValue" and boundary_type != "fixedGradient":
            print "not supported boundary type '" + boundary_type + "' in set_patch_uniform"
            print "supported types are: 'fixedValue' and 'fixedGradient'"
            sys.exit()
            
        self.type = boundary_type 
            
        # check length of boundary value vector
        if len(boundary_values) != len(self.father[0].faces):
            print "define distributed patch with "+str(len(self.father[0].faces))+" faces, only "+str(len(boundary_values))+" defined..."
            sys.exit()
            
        vector_field = False
        if (self.field.type == "vector"):
            vector_field = True
            
        fixed_value = False
        if boundary_type == "fixedValue":
            fixed_value = True
        
        if vector_field and fixed_value:
            self.values = numpy.zeros((len(self.father[0].faces),3))
            for face_i,face_ in enumerate(self.father[0].faces):
                self.values[face_i, :] = boundary_values[face_i]
                face_.bndType = boundary_type
        else:
            self.values = numpy.zeros(len(self.father[0].faces))
            for face_i,face_ in enumerate(self.father[0].faces):
                self.values[face_i] = boundary_values[face_i]
                face_.bndType = boundary_type
                
    def find_face_id_in_patch(self, face_):
        """
        return face id within the patch
        
        :param face_: face object to find
        :type face_:  :class:`pyCFD_mesh.face.face`
        :return:      face id in patch
        :rtype:       int
        """        
        return face_.id - self.father[0].ids[0]