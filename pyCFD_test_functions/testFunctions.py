"""
module of functions for testing mesh data
"""
__author__ = "Bence Somogyi"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Bence Somogyi"
__email__ = "bencesomogyi@ivt.tugraz.at"
__status__ = "Prototype"

import numpy
import pyCFD_VTK_tools.vtkTools as vtkTools

def write_element_face_vector_files(MSHMeshObject):
    for cellI in range(len(MSHMeshObject.cells)):
        cell = MSHMeshObject.cells[cellI]
        fcIter = 0
        for faceI in range(len(cell.faces)):
            face = cell.faces[faceI]
            fileName = '_OUTPUT/e_'+str(cellI)+'_f_'+str(faceI)+'vect.vtu'
            vtkTools.save_vtu_vector(face.C,cell.faceVector[faceI],fileName)
            fcIter += 1
            
def write_internal_face_ffToF_vector_files(mesh_):
    for face_ in mesh_.faces:
        if face_.isBnd:
            continue
        vect_center = numpy.add(-face_.ffToF, face_.C)
        file_name = '_OUTPUT/ffToF_'+str(face_.id)+'vect.vtu'
        vtkTools.save_vtu_vector(vect_center, face_.ffToF, file_name)
        
def write_owner_to_neighbour_vectors(mesh_):
    for face_ in mesh_.faces:
        if face_.isBnd:
            continue
        vect_o_to_n = numpy.add(face_.cells[1].C,-face_.cells[0].C)
        file_name = '_OUTPUT/oToN_'+str(face_.id)+'vect.vtu'
        vtkTools.save_vtu_vector(face_.cells[0].C, vect_o_to_n, file_name)

def write_face_files(MSHMeshObject):
    for faceI in range(len(MSHMeshObject.faces)):
        face = MSHMeshObject.faces[faceI]
        fileName = '_OUTPUT/f_'+str(faceI)+'.vtu'
        coordList = []
        for vertex in face.vertices:
            coordList.append(vertex.get_coords())
        vtkTools.save_vtu_face(coordList,fileName)

def write_patch_faces(MSHMeshObject):
    for patch in MSHMeshObject.patches:
        for faceI in range(len(patch.faces)):
            fileName = '_OUTPUT/'+str(patch.name)+'_f_'+str(faceI)+'.vtu'
            coordList = []
            for vertex in patch.faces[faceI].vertices:
                coordList.append(vertex.get_coords())
            vtkTools.save_vtu_face(coordList,fileName)

def print_mesh_data(meshObject):
    """print minimal info about mesh"""
    print '\n'
    print '===MESH DATA:==='
    print 'number of cells: '+str(len(meshObject.cells))
    print 'number of faces: '+str(len(meshObject.faces))
    print 'number of patches: '+str(len(meshObject.patches))
    print 'patch names: '
    if len(meshObject.patchNames) != 0:
        for patchName in meshObject.patchNames:
            print '\t'+patchName
    print '================'



