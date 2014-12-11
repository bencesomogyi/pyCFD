"""
module of functions for writing data files
"""
__author__ = "Bence Somogyi"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Bence Somogyi"
__email__ = "bencesomogyi@ivt.tugraz.at"
__status__ = "Prototype"

import pyCFD_config.config as Config_
import os
import time
import pyCFD_VTK_tools.vtkTools as vtkTools
import pyCFD_fields.fields as fields
import pyCFD_mesh.sub_mesh as sub_mesh

# "global variables"
__OUTITERDIR__    = Config_.__OUTITERDIR__
__OUTITERDIRREL__ = Config_.__OUTITERDIRREL__
__OUTDIR__        = Config_.__OUTDIR__

def clean_output_dirs():
    """
    delete files from :mod:`pyCFD_config.config`.__OUTDIR__ and
    :mod:`pyCFD_config.config`.__OUTITERDIR__
    """
    for folder in [__OUTITERDIR__, __OUTDIR__]:
        for the_file in os.listdir(folder):
            file_path = os.path.join(folder, the_file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
            except Exception, e:
                print e

def write_patch_files(mesh_):
    """
    write separate .vtu files containing patch faces for all the patches in the
    :mod:`pyCFD_config.config`.__OUTITERDIR__
    
    :param mesh_:   mesh object
    :type mesh_:    :class:`pyCFD_mesh.generic_emsh.GenericMesh`
    """
    t_start = time.time()
    print ''
    print 'Writing patch files...'
    for patch in mesh_.patches:
        fileName = __OUTITERDIR__+str(patch.name)+'.vtu'
        vtkTools.save_vtu_objects(patch.faces,fileName)
        print '   DONE for '+str(patch.name)
    print '   All DONE in '+str(time.time()-t_start)+' s'
    

def write_mesh_file(Mesh):
    """
    write .vtu file containing all the cells in the
    :mod:`pyCFD_config.config`.__OUTITERDIR__.
    
    :param mesh_:   mesh object
    :type mesh_:    :class:`pyCFD_mesh.generic_emsh.GenericMesh`
    """
    t_start = time.time()
    print ''
    print 'Writing mesh file...'
    fileName = __OUTITERDIR__+'mesh.vtu'
    vtkTools.save_vtu_objects(Mesh.cells,fileName)
    print '   DONE in '+str(time.time()-t_start)+' s'

def write_mesh_file_with_fields(field_list, name_=""):
    """
    write .vtu file containing all the cells and field values in the
    :mod:`pyCFD_config.config`.__OUTITERDIR__ directory
    
    :param fieldList:   a list of fields to include in the output. *father* of fields defines the mesh
    :type fieldList:    :class:`pyCFD_fields.fields.VolumeField`
    :param name_:       filename to use
    :type name_:        string
    """
    t_start = time.time()
    print ''
    if name_ == "":
        print 'Writing field file...'
        file_name = __OUTITERDIR__+'fields.vtu'
    else:
        print 'Writing field file '+name_+' ...'
        file_name = __OUTITERDIR__+name_+'.vtu'
    only_one_field = False
    if isinstance(field_list, fields.Field):
        only_one_field = True
    if only_one_field:
        mesh_object = field_list.father[0]
    else:
        mesh_object = field_list[0].father[0]
    field_list_ = []
    if only_one_field:
        field_list_.append(field_list)
    else:
        field_list_.extend(field_list)
    vtkTools.save_vtu_objects(mesh_object.cells, file_name, field_list_)
    print '   DONE in '+str(time.time()-t_start)+' s'
    
def write_patch_files_with_fields(mesh_, volume_field_list):
    """
    write patch meshes with field values into separate .vtu files in the
    :mod:`pyCFD_config.config`.__OUTITERDIR__ directory
    
    :param mesh_:             a mesh object
    :type mesh_:              :class:`pyCFD_mesh.generic_emsh.GenericMesh`
    :param volume_field_list: list of volume fields
    :type volume_field_list:  :class:`pyCFD_fields.fields.VolumeField`
    """
    t_start = time.time()
    print ''
    print 'Writing patch files...'
    only_one_field = False
    if isinstance(volume_field_list, fields.Field):
        only_one_field = True
    field_list_ = []
    if only_one_field:
        field_list_.append(volume_field_list)
    else:
        field_list_.extend(volume_field_list)
    for patch_ in mesh_.patches:
        for field_ in field_list_:
            field_.update_boundary_values()
        fileName = __OUTITERDIR__+str(patch_.name)+'.vtu'
        vtkTools.save_vtu_objects(patch_.faces, fileName, field_list_)
        print '   DONE for '+str(patch_.name)
    print '   All DONE in '+str(time.time()-t_start)+' s'
    
def write_patch_file_with_fields(mesh_, patch_name, volume_field_list):
    """
    write a patch mesh with field values into separate a .vtu file in the
    :mod:`pyCFD_config.config`.__OUTITERDIR__ directory
    
    :param mesh_:             a mesh object
    :type mesh_:              :class:`pyCFD_mesh.generic_emsh.GenericMesh`
    :param patch_name:        name of geometric patch
    :type patch_name:         string
    :param volume_field_list: list of volume fields
    :type volume_field_list:  :class:`pyCFD_fields.fields.VolumeField`
    """
    t_start = time.time()
    print ''
    print 'Writing patch file for '+patch_name
    only_one_field = False
    if isinstance(volume_field_list, fields.Field):
        only_one_field = True
    field_list_ = []
    if only_one_field:
        field_list_.append(volume_field_list)
    else:
        field_list_.extend(volume_field_list)
    patch_ = mesh_.get_patch(patch_name)
    for field_ in field_list_:
        field_.update_boundary_values()
    fileName = __OUTITERDIR__+str(patch_.name)+'.vtu'
    vtkTools.save_vtu_objects(patch_.faces, fileName, field_list_)
    print '   DONE for '+str(patch_.name)+' in '+str(time.time()-t_start)+' s'
    
def write_mesh_internal_faces(mesh_):
    """
    write .vtu file with all the internal faces in the
    :mod:`pyCFD_config.config`.__OUTITERDIR__ directory
    
    :param mesh_:             a mesh object
    :type mesh_:              :class:`pyCFD_mesh.generic_emsh.GenericMesh`
    """
    t_start = time.time()
    print ''
    print 'Writing mesh file with internal faces...'
    fileName = __OUTITERDIR__+'meshIntFaces.vtu'
    internal_faces = []
    for face_ in mesh_.faces:
        if face_.isBnd == False:
            internal_faces.append(face_)
    internal_face_mesh = sub_mesh.SubMesh(None,internal_faces)
    surf_normal = fields.SurfaceVectorField(internal_face_mesh,'Sf')
    for i in range(len(surf_normal.A)):
        surf_normal.A[i] = surf_normal.father[0].faces[i].Sf
    fieldList = []
    fieldList.append(surf_normal)
    vtkTools.save_vtu_objects(internal_face_mesh.faces,fileName,fieldList,mesh_)
    print '   DONE in '+str(time.time()-t_start)+' s'
    
def write_mesh_internal_faces_with_field(surface_field_list):
    """
    write .vtu file with all the internal faces in the
    :mod:`pyCFD_config.config`.__OUTITERDIR__ directory
    
    :param mesh_:             a mesh object
    :type mesh_:              :class:`pyCFD_mesh.generic_emsh.GenericMesh`
    """
    t_start = time.time()
    print ''
    print 'Writing mesh file with internal faces...'
    file_name = __OUTITERDIR__+'meshIntFaces.vtu'
    internal_faces = []
    only_one_field = False
    if isinstance(surface_field_list, fields.Field):
        only_one_field = True
    if only_one_field:
        mesh_object = surface_field_list.father[0]
    else:
        mesh_object = surface_field_list[0].father[0]
    for face_ in mesh_object.faces:
        if face_.isBnd == False:
            internal_faces.append(face_)
    internal_face_mesh = sub_mesh.SubMesh(None,internal_faces)
    field_list = []
    if only_one_field:
        field_list.append(surface_field_list)
    else:
        field_list.extend(surface_field_list)
    vtkTools.save_vtu_objects(internal_face_mesh.faces, file_name, field_list, mesh_object)
    print '   DONE in '+str(time.time()-t_start)+' s'
    
def write_mesh_faces_with_field(surface_field_list, name_=""):
    """
    write .vtu file with all the faces and a list of surface field in the 
    :mod:`pyCFD_config.config`.__OUTITERDIR__ directory
    
    :param mesh_:             a mesh object
    :type mesh_:              :class:`pyCFD_mesh.generic_emsh.GenericMesh`
    """
    t_start = time.time()
    print ''
    if name_ == "":
        file_name = __OUTITERDIR__+'meshAllFaces.vtu'
        print "Writing mesh file 'meshAllFaces.vtu' with all faces..."
    else:
        file_name = __OUTITERDIR__ + name_ + '.vtu'
        print "Writing mesh file '" + name_ + ".vtu' with all faces..."
    only_one_field = False
    if isinstance(surface_field_list, fields.Field):
        only_one_field = True
    if only_one_field:
        mesh_object = surface_field_list.father[0]
    else:
        mesh_object = surface_field_list[0].father[0]
    field_list = []
    if only_one_field:
        field_list.append(surface_field_list)
    else:
        field_list.extend(surface_field_list)
    vtkTools.save_vtu_objects(mesh_object.faces, file_name, field_list, mesh_object)
    print '   DONE in '+str(time.time()-t_start)+' s'
    
def write_mesh_faces(mesh_):
    """
    write .vtu file with all the faces in the 
    :mod:`pyCFD_config.config`.__OUTITERDIR__ directory
    
    :param mesh_:             a mesh object
    :type mesh_:              :class:`pyCFD_mesh.generic_emsh.GenericMesh`
    """
    t_start = time.time()
    print ''
    print 'Writing mesh file with all faces...'
    fileName = __OUTITERDIR__+'meshAllFaces.vtu'
    surf_normal = fields.SurfaceVectorField(mesh_,'Sf')
    for i in range(len(surf_normal.A)):
        surf_normal.A[i] = surf_normal.father[0].faces[i].Sf
    fieldList = []
    fieldList.append(surf_normal)
    vtkTools.save_vtu_objects(mesh_.faces,fileName,fieldList)
    print '   DONE in '+str(time.time()-t_start)+' s'
    
def write_pvd_collection(times_):
    """
    write a .pvd file with links to all the .vtu files including time information
    in the :mod:`pyCFD_config.config`.__OUTDIR__ directory
    
    :param times_: array of timesteps to include
    :type times_:  float
    """
    file_name = 'iterations.pvd'
    vtkTools.save_pvd_collection(times_, __OUTDIR__, __OUTITERDIRREL__, file_name)
    
        
def write_standalone_pvd_collection():
    """
    write a .pvd file with links to all the .vtu files without time information
    in the :mod:`pyCFD_config.config`.__OUTDIR__ directory. All .vtu files in
    the :mod:`pyCFD_config.config`.__OUTITERDIR__ directory will be included.
    """
    file_name = 'iterations.pvd'
    times_ = [i_ for i_ in range(len(os.listdir(__OUTITERDIR__))-1)]
    vtkTools.save_pvd_collection(times_, __OUTDIR__, __OUTITERDIRREL__, file_name, True)
    
def write_csv_file(headers, vectors, file_name_with_path):
    columns_ = len(headers)
    rows_    = len(vectors[0])
    f = open(file_name_with_path, 'w')
    for col_ in range(columns_):
        f.write(headers[col_])
        if col_ != (columns_-1):
            f.write(',')
        else:
            f.write('\n')
    for i_ in range(rows_):
        for col_ in range(columns_):
            f.write(str(vectors[col_][i_]))
            if col_ != (columns_-1):
                f.write(',')
            else:
                f.write('\n')
    f.close()
    
def write_row_to_csv_file(data_row, file_name_with_path, write_or_append = "w"):
    columns_ = len(data_row)
    f = open(file_name_with_path, write_or_append)
    for col_ in range(columns_):
        if write_or_append == 'w':
            append_string = data_row[col_]
        else:
            append_string = str(data_row[col_])
        f.write(append_string)
        if col_ != (columns_-1):
            f.write(',')
        else:
            f.write('\n')
    f.close()