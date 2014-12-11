"""
This module provides functions to write ascii VTK or VTU files that can be
later loaded for postprocessing in Paraview.
"""
__author__ = "Bence Somogyi"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Bence Somogyi"
__email__ = "bencesomogyi@ivt.tugraz.at"
__status__ = "Prototype"

import os
import sys
import pyCFD_mesh.face as meshFace
import pyCFD_mesh.cell as meshCell
import pyCFD_mesh.vertex as meshVertex

def save_vtu_vector(pointCoords,vectCoords,fileNameWithPath):
    """
    function to write VTU file of a vector at a point
    
    :param pointCoords:      list of vertex coordinates
    :type pointCoords:       float
    :param vectCoords:       list of vector coordinates
    :type vectCoords:        float
    :param fileNameWithPath: location and file name to save
    :type fileNameWithPath:  string
    """
    f = open(fileNameWithPath, 'w')
    f.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
    f.write('\t<UnstructuredGrid>\n')
    f.write('\t\t<Piece NumberOfPoints="1" NumberOfCells="1">\n')
    f.write('\t\t\t<PointData Vectors="vector">\n')
    f.write('\t\t\t\t<DataArray type="Float32" Name="vector" NumberOfComponents="3" format="ascii">\n')
    f.write('\t\t\t\t\t')
    for i in range(3):
        f.write(' ')
        f.write(str(vectCoords[i]))
    f.write('\n')
    f.write('\t\t\t\t</DataArray>\n')
    f.write('\t\t\t</PointData>\n')
    f.write('\t\t\t<Points>\n')
    f.write('\t\t\t\t<DataArray type="Float32" NumberOfComponents="3" format="ascii">\n')
    f.write('\t\t\t\t\t')
    for i in range(3):
        f.write(' ')
        f.write(str(pointCoords[i]))
    f.write('\n')
    f.write('\t\t\t\t</DataArray>\n')
    f.write('\t\t\t</Points>\n')
    f.write('\t\t\t<Cells>\n')
    f.write('\t\t\t\t<DataArray type="Int32" Name="connectivity" format="ascii">\n')
    f.write('\t\t\t\t\t0\n')
    f.write('\t\t\t\t</DataArray>\n')
    f.write('\t\t\t\t<DataArray type="Int32" Name="offsets" format="ascii">\n')
    f.write('\t\t\t\t\t1\n')
    f.write('\t\t\t\t</DataArray>\n')
    f.write('\t\t\t\t<DataArray type="Int32" Name="types" format="ascii">\n')
    f.write('\t\t\t\t\t1\n')
    f.write('\t\t\t\t</DataArray>\n')
    f.write('\t\t\t</Cells>\n')
    f.write('\t\t</Piece>\n')
    f.write('\t</UnstructuredGrid>\n')
    f.write('</VTKFile>')
    f.close()

def save_vtu_face(pointCoords,fileNameWithPath):
    """
    function to write VTU file of a face
    
    :param pointCoords:      list of vertex coordinates
    :type volume_field:      float
    :param fileNameWithPath: location and file name to save
    :type fileNameWithPath:  string
    """
    if len(pointCoords) == 4:
        # vtk quadrangle
        vtkType = 9
    elif len(pointCoords) == 3:
        # vtk triangle
        vtkType = 5
    else:
        sys.exit("saving this type of VTU face is not supported!")

    f = open(fileNameWithPath, 'w')
    f.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
    f.write('\t<UnstructuredGrid>\n')
    f.write('\t\t<Piece NumberOfPoints="'+str(len(pointCoords))+'" NumberOfCells="1">\n')
    f.write('\t\t\t<Points>\n')
    f.write('\t\t\t\t<DataArray type="Float32" NumberOfComponents="3" format="ascii">\n')
    for i in range(len(pointCoords)):
        f.write('\t\t\t\t\t')
        for j in range(3):
            f.write(str(pointCoords[i][j]))
            f.write(' ')
        f.write('\n')
    f.write('\t\t\t\t</DataArray>\n')
    f.write('\t\t\t</Points>\n')
    f.write('\t\t\t<Cells>\n')
    f.write('\t\t\t\t<DataArray type="Int32" Name="connectivity" format="ascii">\n')
    f.write('\t\t\t\t\t')
    for i in range(len(pointCoords)):
        f.write(str(i)+" ")
    f.write('\n')
    f.write('\t\t\t\t</DataArray>\n')
    f.write('\t\t\t\t<DataArray type="Int32" Name="offsets" format="ascii">\n')
    f.write('\t\t\t\t\t'+str(len(pointCoords))+'\n')
    f.write('\t\t\t\t</DataArray>\n')
    f.write('\t\t\t\t<DataArray type="Int32" Name="types" format="ascii">\n')
    f.write('\t\t\t\t\t'+str(vtkType)+'\n')
    f.write('\t\t\t\t</DataArray>\n')
    f.write('\t\t\t</Cells>\n')
    f.write('\t\t</Piece>\n')
    f.write('\t</UnstructuredGrid>\n')
    f.write('</VTKFile>')
    f.close()

def save_vtu_objects(objList,fileNameWithPath,fieldList=None,fatherMesh=None,writeVector=False):
    """
    function to write VTK file with objects (faces/cells) given in objList
    
    :param objList:          list of mesh objects
    :type objList:           :class:`pyCFD_mesh.mesh_object.MeshObject`
    :param fileNameWithPath: location and file name to save
    :type fileNameWithPath:  string
    :param fieldList:        default None, list of volume fields to save
    :type fieldList:         :class:`pyCFD_fields.fields.VolumeField`
    :param fatherMesh:       default: None, mesh objects' owner mesh
    :type fatherMesh:        :class:`pyCFD_mesh.generic_mesh.GenericMesh`
    :param writeVector:      default: False, if vector fields should be written
    :type writeVector:       bool
    """
    vertices    = []
    newVertices = []
    vtkTypes    = []
    offsets     = []
    objIds      = []
    scalList    = []
    vectList    = []

    def fill_obj_data_lists():
        # faces or tetrahedron
        vertex_number = len(obj.vertices)
        if faceMode or vertex_number < 6:
            newVertCurrent = []
            for vertex_ in obj.vertices:
                newVertCurrent.append(vertex_.id)
            newVertices.append(newVertCurrent)
        else:
            newVertCurrent = []
            if vertex_number == 6: # prism
                for face_ in obj.faces:
                    if len(face_.vertices) == 3:
                        face_1 = face_
            else: # hexahedron
                face_1 = obj.faces[0]
            newVertCurrent.extend(face_1.vertices)
            vertex_1 = face_1.vertices[0]
            # find face2 vertices
            face_2_vertices = []
            for vert2 in obj.vertices:
                if vert2 not in face_1.vertices:
                    face_2_vertices.append(vert2)
            # find face2
            sorted_vert_face_2 = meshVertex.get_list_of_ids(face_2_vertices)
            sorted_vert_face_2.sort()
            for face_ in obj.faces:
                sorted_vert_face_ = face_.get_vertex_ids()
                sorted_vert_face_.sort()
                if face_ == face_1:
                    continue
                if sorted_vert_face_ == sorted_vert_face_2:
                    face_2 = face_
                    break
            # find connecting face
            connecting_face = []
            for face_ in vertex_1.faces:
                if face_ == face_1:
                    continue
                cell_ids = face_.get_cell_ids()
                if obj.id not in cell_ids:
                    continue
                connecting_face.append(face_)
            # find common vertices of conn_face1/face_2 and conn_face2/face_2
            common_vertices_0 = []
            common_vertices_1 = []
            for vertex_ in connecting_face[0].vertices:
                if vertex_ in face_2.vertices:
                    common_vertices_0.append(vertex_)
            for vertex_ in connecting_face[1].vertices:
                if vertex_ in face_2.vertices:
                    common_vertices_1.append(vertex_)
            # find neighbours of vertex_1
            face_2_vert_2 = []
            face_2_vert_2.extend(face_2.vertices)
            neigbours_to_vertex_1_0 = connecting_face[0].find_neighbour_vertex_on_face(vertex_1)
            for vertex_ in neigbours_to_vertex_1_0:
                if vertex_ in common_vertices_0:
                    newVertCurrent.append(vertex_)
                    break
            common_vertices_0.remove(newVertCurrent[-1])
            common_vertices_1.remove(newVertCurrent[-1])
            face_2_vert_2.remove(common_vertices_0[0])
            face_2_vert_2.remove(common_vertices_1[0])
            face_2_vert_2.remove(newVertCurrent[-1])
            if vertex_number == 6:
                if newVertCurrent[1] in connecting_face[0].vertices:
                    newVertCurrent.append(common_vertices_0[0])
                    newVertCurrent.append(common_vertices_1[0])
                else:
                    newVertCurrent.append(common_vertices_1[0])
                    newVertCurrent.append(common_vertices_0[0])
            else:                    
                if newVertCurrent[1] in connecting_face[0].vertices:
                    newVertCurrent.append(common_vertices_0[0])
                    newVertCurrent.append(face_2_vert_2[0])
                    newVertCurrent.append(common_vertices_1[0])
                else:
                    newVertCurrent.append(common_vertices_1[0])
                    newVertCurrent.append(face_2_vert_2[0])
                    newVertCurrent.append(common_vertices_0[0])
            newVertices.append(meshVertex.get_list_of_ids(newVertCurrent))
        # create list with original obj ids from the mesh
        if len(fullObjList) == len(objList):
            objInd = obj.id
        else:
            objInd = fullObjList.index(obj)
        objIds.append(objInd)
        pass

    def add_fields(fileStream,_fieldList):
        if _fieldList == None:
            pass
        else:
            for field in _fieldList:
                if field.type == 'scalar':
                    scalList.append(field)
                if field.type == 'vector':
                    vectList.append(field)
                    
            # CellData
            fileStream.write('\t\t\t<CellData')
            if len(scalList) != 0:
                fileStream.write(' Scalars="')
                for scalField in scalList:
                    fileStream.write(' '+str(scalField.name))
                fileStream.write('"')
            if len(vectList) != 0:
                fileStream.write(' Vectors="')
                for vectField in vectList:
                    fileStream.write(' '+str(vectField.name))
                fileStream.write('"')
            fileStream.write('>\n')
            for scalField in scalList:
                fileStream.write('\t\t\t\t<DataArray type="Float32" Name="'+str(scalField.name)+'" format="ascii">\n')
                if cellMode:
                    sourceField = scalField.V
                elif faceMode:
                    sourceField = scalField.A
                for objI in range(len(objList)):
                    fileStream.write('\t\t\t\t\t')
                    fileStream.write(str(sourceField[objIds[objI]]))
                    fileStream.write('\n')
                f.write('\t\t\t\t</DataArray>\n')
            for vectField in vectList:
                fileStream.write('\t\t\t\t<DataArray type="Float32" NumberOfComponents="3" Name="'+str(vectField.name)+'" format="ascii">\n')
                if cellMode:
                    sourceField = vectField.V
                elif faceMode:
                    sourceField = vectField.A
                for objI in range(len(objList)):
                    fileStream.write('\t\t\t\t\t')
                    fileStream.write(str(sourceField[objIds[objI]][0]))
                    fileStream.write(' ')
                    fileStream.write(str(sourceField[objIds[objI]][1]))
                    fileStream.write(' ')
                    fileStream.write(str(sourceField[objIds[objI]][2]))
                    fileStream.write('\n')
                fileStream.write('\t\t\t\t</DataArray>\n')
            fileStream.write('\t\t\t</CellData>\n')
            # PointData in separated file

    # set mode: cell/face
    if isinstance(objList[0],meshFace.Face):
        faceMode = True
    if isinstance(objList[0],meshCell.Cell):       
        faceMode = False
    cellMode = not faceMode

    # check consistency of objects
    for obj in objList[1:]:
        if faceMode:
            if isinstance(obj,meshFace.Face) is False:
                sys.exit("inconsistent object list!")
        if cellMode:
            if isinstance(obj,meshCell.Cell) is False:
                sys.exit("inconsistent object list!")

    # list of original data
    if fatherMesh == None:
        father = objList[0].father[0]
    else:
        father = fatherMesh
        
    vertices = father.vertices
        
    if cellMode:
        fullObjList = father.cells
    elif faceMode:
        fullObjList = father.faces

    isFirstObj = True
    if faceMode:
        for objI,obj in enumerate(objList):
            # find face types and offsets
            if len(obj.vertices) == 4:
                # vtk quadrangle
                vtkTypes.append(9)
                if isFirstObj:
                    offsets.append(4)
                else:
                    offsets.append(offsets[-1]+4)
            elif len(obj.vertices) == 3:
                # vtk triangle
                vtkTypes.append(5)
                if isFirstObj:
                    offsets.append(3)
                else:
                    offsets.append(offsets[-1]+3)
            else:
                print "face id: "+str(obj.id)+" with "+str(len(obj.vertices))+" faces"
                sys.exit("saving this type of VTU face is not supported!")
            fill_obj_data_lists()
            isFirstObj = False
    elif cellMode:
        for objI,obj in enumerate(objList):
            # find cell types and offsets
            if len(obj.vertices) == 4:
                # vtk tetra
                vtkTypes.append(10)
                if isFirstObj:
                    offsets.append(4)
                else:
                    offsets.append(offsets[-1]+4)
            elif len(obj.vertices) == 6:
                # vtk wedge
                vtkTypes.append(13)
                if isFirstObj:
                    offsets.append(6)
                else:
                    offsets.append(offsets[-1]+6)
            elif len(obj.vertices) == 8:
                # vtk hexahedron
                vtkTypes.append(12)
                if isFirstObj:
                    offsets.append(8)
                else:
                    offsets.append(offsets[-1]+8)
            else:
                sys.exit("saving this type of VTU cell is not supported!")
            fill_obj_data_lists()
            isFirstObj = False

    f = open(fileNameWithPath, 'w')
    f.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
    f.write('\t<UnstructuredGrid>\n')
    f.write('\t\t<Piece NumberOfPoints="'+str(len(vertices))+'" NumberOfCells="'+str(len(objList))+'">\n')
    add_fields(f,fieldList)
    f.write('\t\t\t<Points>\n')
    f.write('\t\t\t\t<DataArray type="Float32" NumberOfComponents="3" format="ascii">\n')
    for vertex in vertices:
        f.write('\t\t\t\t\t')
        f.write(str(vertex.X))
        f.write(' ')
        f.write(str(vertex.Y))
        f.write(' ')
        f.write(str(vertex.Z))
        f.write('\n')
    f.write('\t\t\t\t</DataArray>\n')
    f.write('\t\t\t</Points>\n')
    f.write('\t\t\t<Cells>\n')
    f.write('\t\t\t\t<DataArray type="Int32" Name="connectivity" format="ascii">\n')
    for objVertexList in newVertices:
        f.write('\t\t\t\t\t')
        for objVertex in objVertexList:
            f.write(str(objVertex)+" ")
        f.write('\n')
    f.write('\t\t\t\t</DataArray>\n')
    f.write('\t\t\t\t<DataArray type="Int32" Name="offsets" format="ascii">\n')
    for objOffset in offsets:
        f.write('\t\t\t\t\t'+str(objOffset)+'\n')
    f.write('\t\t\t\t</DataArray>\n')
    f.write('\t\t\t\t<DataArray type="Int32" Name="types" format="ascii">\n')
    for objType in vtkTypes:
        f.write('\t\t\t\t\t'+str(objType)+'\n')
    f.write('\t\t\t\t</DataArray>\n')
    f.write('\t\t\t</Cells>\n')
    f.write('\t\t</Piece>\n')
    f.write('\t</UnstructuredGrid>\n')
    f.write('</VTKFile>')
    f.close()

    if writeVector:
        save_vtu_vector_field(objList,vectList,fileNameWithPath)

def save_vtu_vector_field(objList,vectList,fileNameWithPath):
    """
    function to write VTU with vector fields
    
    :param objList:          list of mesh objects
    :type objList:           :class:`pyCFD_mesh.mesh_object.MeshObject`
    :param vectList:         list of volume vector fields to save
    :type vectList:          :class:`pyCFD_fields.fields.VectorField`
    :param fileNameWithPath: location and file name to save
    :type fileNameWithPath:  string
    """
    if isinstance(objList[0],meshFace.Face):
        faceMode = True
    if isinstance(objList[0],meshCell.Cell):
        faceMode = False
    cellMode = not faceMode

    vertices = [None]*len(objList)
    objIds   = [None]*len(objList)

    for i,obj in enumerate(objList):
        vertices[i] = meshVertex.Vertex(obj.C[0],obj.C[1],obj.C[2])
        if cellMode:
            objInd = obj.father[0].cells.index(obj)
        elif faceMode:
            objInd = obj.father[0].faces.index(obj)
        objIds[i] = objInd

    newFileNameWithPath = os.path.splitext(fileNameWithPath)[0] + '_vector' + os.path.splitext(fileNameWithPath)[1]

    def add_fields(fileStream,_vectList):
        # PointData
        fileStream.write('\t\t\t<PointData')
        if len(vectList) != 0:
            fileStream.write(' Vectors="')
            for vectField in vectList:
                fileStream.write(' '+str(vectField.name))
            fileStream.write('"')
        fileStream.write('>\n')
        for vectField in vectList:
            fileStream.write('\t\t\t\t<DataArray type="Float32" NumberOfComponents="3" Name="'+str(vectField.name)+'" format="ascii">\n')
            if cellMode:
                sourceField = vectField.V
            elif faceMode:
                sourceField = vectField.A
            for objI in range(len(objList)):
                fileStream.write('\t\t\t\t\t')
                fileStream.write(str(sourceField[objIds[objI]][0]))
                fileStream.write(' ')
                fileStream.write(str(sourceField[objIds[objI]][1]))
                fileStream.write(' ')
                fileStream.write(str(sourceField[objIds[objI]][2]))
                fileStream.write('\n')
            fileStream.write('\t\t\t\t</DataArray>\n')
        fileStream.write('\t\t\t</PointData>\n')

    f = open(newFileNameWithPath, 'w')
    f.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
    f.write('\t<UnstructuredGrid>\n')
    f.write('\t\t<Piece NumberOfPoints="'+str(len(vertices))+'" NumberOfCells="0">\n')
    add_fields(f,vectList)
    f.write('\t\t\t<Points>\n')
    f.write('\t\t\t\t<DataArray type="Float32" NumberOfComponents="3" format="ascii">\n')
    for vertex in vertices:
        f.write('\t\t\t\t\t')
        f.write(str(vertex.X))
        f.write(' ')
        f.write(str(vertex.Y))
        f.write(' ')
        f.write(str(vertex.Z))
        f.write('\n')
    f.write('\t\t\t\t</DataArray>\n')
    f.write('\t\t\t</Points>\n')
    f.write('\t\t\t<Cells>\n')
    f.write('\t\t\t\t<DataArray type="Int32" Name="connectivity" format="ascii">\n')
    f.write('\t\t\t\t</DataArray>\n')
    f.write('\t\t\t\t<DataArray type="Int32" Name="offsets" format="ascii">\n')
    f.write('\t\t\t\t</DataArray>\n')
    f.write('\t\t\t\t<DataArray type="Int32" Name="types" format="ascii">\n')
    f.write('\t\t\t\t</DataArray>\n')
    f.write('\t\t\t</Cells>\n')
    f.write('\t\t</Piece>\n')
    f.write('\t</UnstructuredGrid>\n')
    f.write('</VTKFile>')
    f.close()

def save_pvd_collection(times_, save_path, data_path, file_name, standalone=False):
    """
    function to PVD with time information and references to the saved VTU files
    
    :param times_:           list of timesteps to save
    :type times_:            list of float
    :param save_path:        location to save PVD
    :type save_path:         string
    :param data_path:        location where VTU files are stored
    :type data_path:         string
    :param standalone:       default: False, if True time information is skipped,
                             all VTU files from the target directory are saved
                             in the PVD
    :type vectList:          bool
    """
    f = open(save_path+"/"+file_name, 'w')
    f.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n')
    f.write('\t<Collection>\n')
    if standalone == False:
        if os.path.isfile(save_path+data_path+"0.vtu"):
            f.write('\t\t<DataSet timestep="0.0" group="" part="0" file="'+data_path+'0.vtu"/>\n')
    for i_,time_ in enumerate(times_):
        f.write('\t\t<DataSet timestep="'+str(time_)+'" group="" part="0" file="'+data_path+str(i_+1)+'.vtu"/>\n')
    f.write('\t</Collection>\n')
    f.write('</VTKFile>')
    f.close()

