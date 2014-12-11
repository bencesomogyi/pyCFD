"""
module for geometric calculations as face areas, cell volumes and centroid
coordinates
"""
__author__ = "Bence Somogyi"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Bence Somogyi"
__email__ = "bencesomogyi@ivt.tugraz.at"
__status__ = "Prototype"

import numpy
import math
import sys
import pyCFD_operators.vector_operations

if sys.platform == 'win32':
    import pyCFD_geometric_tools.cython_boost_win32.cy_geometric_tools as cy_geometric_tools
elif sys.platform == 'linux2':
    import pyCFD_geometric_tools.cython_boost_linux2.cy_geometric_tools as cy_geometric_tools
else:
    sys.exit("unknown platform " + sys.platform + " found in geomTools.py, stopping...")

def triangleCentroid(ndCoordList):
    """
    calculate coordinates of the centroid of a 3D triangle from list of
    coordinates of its vertices
    
    :param ndCoordList: list of arrays with node coordinates
    :type ndCoordList:  numpy.array
    :return:            array with centroid coordinate
    :type:              numpy.array
    """
    [pointACoords, pointBCoords, pointCCoords] = ndCoordList
    centroid = numpy.add(numpy.add(pointACoords,pointBCoords),pointCCoords) / 3.0
    return centroid

def triangle_centroid_vert(vertList):
    """
    calculate coordinates of the centroid of a 3D triangle from a list of
    vertex objects
    
    :param vertList: list of vertex objects
    :type vertList:  :class:`pyCFD_mesh.vertex.Vertex`
    :return:         array with centroid coordinate
    :type:           numpy.array
    """
    pointACoords = vertList[0].coords
    pointBCoords = vertList[1].coords
    pointCCoords = vertList[2].coords
    centroid = numpy.add(numpy.add(pointACoords,pointBCoords),pointCCoords) / 3.0
    return centroid

def polygon_centroid_and_area(vertex_list):
    """
    calculate the centroid and area of a polygon from the list of its vertices
    
    **Calculation steps:**
        
    * polygon is subdivided into triangles using the geometric centroid of the polygon
    
    * centers of the sub-triangles and their area is calculated
    
    * polygon centroid is calculated as surface weighted average of the sub triangles
    
    * polygon area is calculated as the sum of sub triangle areas
    
    :param vertex_list: list of vertex objects defining the face
    :type vertex_list:  :class:`pyCFD_mesh.vertex.Vertex`
    :return:            dictionary [weighted centroid, area]
    :rtype:             float
    """
    number_of_vertices = len(vertex_list)
    if number_of_vertices < 3:
        sys.exit("invalid face element with low number of vertices (< 3)!")
    # get sub triangle coordinates
    vertex_coords = numpy.array([vertex_list[i_vertex].coords for i_vertex in xrange(number_of_vertices)])
    centroid = numpy.array([vertex_coords.sum(0) / float(number_of_vertices),] * number_of_vertices)
    coords_a = numpy.array([vertex_coords[i_vertex] for i_vertex in xrange(number_of_vertices)])
    coords_b = numpy.roll(coords_a,-1,0)
    # calculate sub tirangle centers and areas
    triangle_centers = numpy.add(numpy.add(centroid, coords_a), coords_b) /3.
    # area_vector = numpy.cross(numpy.add(centroid, -coords_a), numpy.add(centroid,-coords_b)) * 0.5
    area_vector = cy_geometric_tools.cy_area_vector(centroid, coords_a, coords_b)
    # triangle_areas = numpy.array([numpy.linalg.norm(area_vector[i]) for i in xrange(number_of_vertices)])
    triangle_areas = cy_geometric_tools.cy_triangle_areas(area_vector)
    # calculate poligon weighted centroid and area    
    weighted_centroid = numpy.array([triangle_centers[i] * triangle_areas[i] for i in xrange(number_of_vertices)])
    area_sum = triangle_areas.sum()
    weighted_centroid_sum = weighted_centroid.sum(0) / area_sum
    result = [weighted_centroid_sum, area_sum]
    return result        

def quadrangleCentroid(ndCoordList):
    """
    calculate coordinates of the centroid of a 3D quadrangle from list of
    coordinates of its vertices
    
    :param ndCoordList: list of arrays with node coordinates
    :type ndCoordList:  numpy.array
    :return:            array with centroid coordinate
    :type:              numpy.array
    """
    [pointACoords, pointBCoords, pointCCoords, pointDCoords] = ndCoordList
    centroid = numpy.add(numpy.add(pointACoords,pointBCoords),numpy.add(pointCCoords,pointDCoords)) / 4.0
    return centroid

def quadrangle_centroid_vert(vertList):
    """
    calculate coordinates of the centroid of a 3D quadrangle from a list of
    vertex objects
    
    :param vertList: list of vertex objects
    :type vertList:  :class:`pyCFD_mesh.vertex.Vertex`
    :return:         array with centroid coordinate
    :type:           numpy.array
    """
    pointACoords = vertList[0].coords
    pointBCoords = vertList[1].coords
    pointCCoords = vertList[2].coords
    pointDCoords = vertList[3].coords
    centroid = numpy.add(numpy.add(pointACoords,pointBCoords),numpy.add(pointCCoords,pointDCoords)) / 4.0
    return centroid

def tetrahedronCentroid(ndCoordList):
    """
    calculate coordinates of the centroid of a tetrahedron from list of
    coordinates of its vertices
    
    :param ndCoordList: list of arrays with node coordinates
    :type ndCoordList:  numpy.array
    :return:            array with centroid coordinate
    :type:              numpy.array
    """
    [pointACoords,pointBCoords,pointCCoords,pointDCoords] = ndCoordList
    centroid = numpy.add(numpy.add(pointACoords,pointBCoords),numpy.add(pointCCoords,pointDCoords)) / 4.0
    return centroid

def tetrahedron_centroid_vert(vertList):
    """
    calculate coordinates of the centroid of a tetrahedron from a list of
    vertex objects
    
    :param vertList: list of vertex objects
    :type vertList:  :class:`pyCFD_mesh.vertex.Vertex`
    :return:         array with centroid coordinate
    :type:           numpy.array
    """
    pointACoords = vertList[0].coords
    pointBCoords = vertList[1].coords
    pointCCoords = vertList[2].coords
    pointDCoords = vertList[3].coords
    centroid = numpy.add(numpy.add(pointACoords,pointBCoords),numpy.add(pointCCoords,pointDCoords)) / 4.0
    return centroid

def hexahedronCentroid(ndCoordList):
    """
    calculate coordinates of the centroid of a hexahedron from list of
    coordinates of its vertices
    
    :param ndCoordList: list of arrays with node coordinates
    :type ndCoordList:  numpy.array
    :return:            array with centroid coordinate
    :type:              numpy.array
    """
    [pointACoords,pointBCoords,pointCCoords,pointDCoords,pointECoords,pointFCoords,pointGCoords,pointHCoords] = ndCoordList
    centroid = numpy.add(numpy.add(numpy.add(pointACoords,pointBCoords),numpy.add(pointCCoords,pointDCoords)),
    numpy.add(numpy.add(pointECoords,pointFCoords),numpy.add(pointGCoords,pointHCoords))) / 8.0
    return centroid

def prism_centroid_vert(vertList):
    """
    calculate coordinates of the centroid of a prism from a list of vertex
    objects
    
    :param vertList: list of vertex objects
    :type vertList:  :class:`pyCFD_mesh.vertex.Vertex`
    :return:         array with centroid coordinate
    :type:           numpy.array
    """
    pointACoords = vertList[0].coords
    pointBCoords = vertList[1].coords
    pointCCoords = vertList[2].coords
    pointDCoords = vertList[3].coords
    pointECoords = vertList[4].coords
    pointFCoords = vertList[5].coords
    centroid = numpy.add(numpy.add(numpy.add(pointACoords,pointBCoords),numpy.add(pointCCoords,pointDCoords)),
    numpy.add(pointECoords,pointFCoords)) / 6.0
    return centroid

def hexahedron_centroid_vert(vertList):
    """
    calculate coordinates of the centroid of a hexahedron from a list of
    vertex objects
    
    :param vertList: list of vertex objects
    :type vertList:  :class:`pyCFD_mesh.vertex.Vertex`
    :return:         array with centroid coordinate
    :type:           numpy.array
    """
    pointACoords = vertList[0].coords
    pointBCoords = vertList[1].coords
    pointCCoords = vertList[2].coords
    pointDCoords = vertList[3].coords
    pointECoords = vertList[4].coords
    pointFCoords = vertList[5].coords
    pointGCoords = vertList[6].coords
    pointHCoords = vertList[7].coords
    centroid = numpy.add(numpy.add(numpy.add(pointACoords,pointBCoords),numpy.add(pointCCoords,pointDCoords)),
    numpy.add(numpy.add(pointECoords,pointFCoords),numpy.add(pointGCoords,pointHCoords))) / 8.0
    return centroid

def triangleAreaVect(ndCoordList, centroidOfElement):
    """
    calculate the face normal area vector of a triangle with respect to the 
    centroid of an element: face normal points outwards
    
    :param ndCoordList:       list of arrays with node coordinates
    :type ndCoordList:        numpy.array
    :param centroidOfElement: array with cell centroid coordinate
    :type centroidOfElement:  numpy.array
    :return:                  array of face area vector
    :type:                    numpy.array
    """
    [pointACoords, pointBCoords, pointCCoords] = ndCoordList
    areaVect = numpy.cross(numpy.add(pointBCoords,-pointACoords), numpy.add(pointCCoords,-pointACoords)) * 0.5
    centroidFace = triangleCentroid(ndCoordList)
    vectElmFace = numpy.add(centroidFace,-centroidOfElement)
    angle = math.acos(numpy.dot(areaVect,vectElmFace) / numpy.linalg.norm(areaVect) / numpy.linalg.norm(vectElmFace)) * 180 / math.pi
    if angle > 90:
        areaVect *= -1
    return areaVect

def face_area_vect_vert(vert_list, centroid_of_element):
    """
    calculate the face normal area vector of a face with respect to the 
    centroid of an element: face normal points outwards
    
    Depending of the number of faces :func:`triangle_area_vect_vert` or 
    func:`quadrangle_area_vect_vert` is called. Other faces are not yet
    supported.
    
    :param vertList:          list of vertex objects
    :type vertList:           :class:`pyCFD_mesh.vertex.Vertex`
    :param centroidOfElement: array with cell centroid coordinate
    :type centroidOfElement:  numpy.array
    :return:                  array of face area vector
    :type:                    numpy.array
    """
    if len(vert_list) == 3:
        return triangle_area_vect_vert(vert_list, centroid_of_element)
    elif len(vert_list) == 4:
        return quadrangle_area_vect_vert(vert_list, centroid_of_element)
    else:
        sys.exit("only triangles and quadrangles are supported. Stopping in geomTools.face_area_vect_vert")

def triangle_area_vect_vert(vertList, centroidOfElement):
    """
    calculate the face normal area vector of a triangle with respect to the 
    centroid of an element: face normal points outwards
    
    :param vertList:          list of vertex objects
    :type vertList:           :class:`pyCFD_mesh.vertex.Vertex`
    :param centroidOfElement: array with cell centroid coordinate
    :type centroidOfElement:  numpy.array
    :return:                  array of face area vector
    :type:                    numpy.array
    """
    pointACoords = vertList[0].coords
    pointBCoords = vertList[1].coords
    pointCCoords = vertList[2].coords
    areaVect = numpy.cross(numpy.add(pointBCoords,-pointACoords), numpy.add(pointCCoords,-pointACoords)) * 0.5
    centroidFace = triangle_centroid_vert(vertList)
    vectElmFace = numpy.add(centroidFace,-centroidOfElement)
    cos_value = numpy.dot(areaVect,vectElmFace) / numpy.linalg.norm(areaVect) / numpy.linalg.norm(vectElmFace)
    if numpy.allclose(cos_value,1.):
        cos_value = 1.
    elif numpy.allclose(cos_value,-1.):
        cos_value = -1.
    angle = math.acos(cos_value) * 180 / math.pi
    if angle > 90:
        areaVect *= -1
    return areaVect
    
def triangle_area_vect_face(face_, centroid_of_element):
    """
    calculate the face normal area vector of a triangle with respect to the 
    centroid of an element: face normal points outwards
    
    :param face_:             a face object
    :type face_:              :class:`pyCFD_mesh.face.Face`
    :param centroidOfElement: array with cell centroid coordinate
    :type centroidOfElement:  numpy.array
    :return:                  array of face area vector
    :type:                    numpy.array
    """
    vertex_list = []
    vertex_list.append(face_.vertices[0])
    vertex_list.append(face_.vertices[1])
    vertex_list.append(face_.vertices[2])
    return triangle_area_vect_vert(vertex_list, centroid_of_element)

def quadrangleAreaVect(ndCoordList, centroidOfElement):
    """
    calculate the face normal area vector of a quadrangle with respect to the 
    centroid of an element: face normal points outwards
    
    :param ndCoordList:       list of arrays with node coordinates
    :type ndCoordList:        numpy.array
    :param centroidOfElement: array with cell centroid coordinate
    :type centroidOfElement:  numpy.array
    :return:                  array of face area vector
    :type:                    numpy.array
    """
    [pointACoords, pointBCoords, pointCCoords, pointDCoords] = ndCoordList
    areaVect = numpy.cross(numpy.add(pointCCoords,-pointACoords), numpy.add(pointBCoords,-pointDCoords)) * 0.5
    centroidFace = quadrangleCentroid(ndCoordList)
    vectElmFace = numpy.add(centroidFace,-centroidOfElement)
    angle = math.acos(numpy.dot(areaVect,vectElmFace) / numpy.linalg.norm(areaVect) / numpy.linalg.norm(vectElmFace)) * 180 / math.pi
    if angle > 90:
        areaVect *= -1
    return areaVect

def quadrangle_area_vect_vert(vertList, centroidOfElement):
    """
    calculate the face normal area vector of a quadrangle with respect to the 
    centroid of an element: face normal points outwards
    
    :param vertList:          list of vertex objects
    :type vertList:           :class:`pyCFD_mesh.vertex.Vertex`
    :param centroidOfElement: array with cell centroid coordinate
    :type centroidOfElement:  numpy.array
    :return:                  array of face area vector
    :type:                    numpy.array
    """
    pointACoords = vertList[0].coords
    pointBCoords = vertList[1].coords
    pointCCoords = vertList[2].coords
    pointDCoords = vertList[3].coords
    areaVect = pyCFD_operators.vector_operations.VecCross(numpy.add(pointCCoords,-pointACoords), numpy.add(pointBCoords,-pointDCoords)) * 0.5
    centroidFace = numpy.add(numpy.add(pointACoords,pointBCoords),numpy.add(pointCCoords,pointDCoords)) / 4.0
    vectElmFace = numpy.add(centroidFace,-centroidOfElement)
#    cos_value = numpy.dot(areaVect,vectElmFace) / numpy.linalg.norm(areaVect) / numpy.linalg.norm(vectElmFace)
    cos_value = cy_geometric_tools.cy_calc_cos(areaVect, vectElmFace)
#    if numpy.allclose(abs(cos_value),1.):
    if cy_geometric_tools.cy_value_close(abs(cos_value),1.):
        if cos_value > 0.:
            cos_value = 1.
        else:
            cos_value = -1.
    angle = math.acos(cos_value) * 180 / math.pi
    if angle > 90:
        areaVect *= -1
    return areaVect
    
def quadrangle_area_vect_face(face_, centroid_of_element):
    """
    calculate the face normal area vector of a quadrangle with respect to the 
    centroid of an element: face normal points outwards
    
    :param face_:             a face object
    :type face_:              :class:`pyCFD_mesh.face.Face`
    :param centroidOfElement: array with cell centroid coordinate
    :type centroidOfElement:  numpy.array
    :return:                  array of face area vector
    :type:                    numpy.array
    """
    vertex_list = []
    vertex_list.append(face_.vertices[0])
    vertex_list.append(face_.vertices[1])
    vertex_list.append(face_.vertices[2])
    vertex_list.append(face_.vertices[3])
    return quadrangle_area_vect_vert(vertex_list, centroid_of_element)

def tetrahedronVolume(ndCoordList):
    """
    the volume of a tetrahedron from the coordinates of its vertices
    
    :param ndCoordList:       list of arrays with node coordinates
    :type ndCoordList:        numpy.array
    :return:                  volume of tetrahedron
    :type:                    float
    """
    [pointACoords,pointBCoords,pointCCoords,pointDCoords] = ndCoordList
    vectDA = numpy.add(pointACoords,-pointDCoords)
    vectDB = numpy.add(pointBCoords,-pointDCoords)
    vectDC = numpy.add(pointCCoords,-pointDCoords)
    volume = math.fabs(numpy.dot(vectDA,numpy.cross(vectDB,vectDC)) / 6)
    return volume
    
def hexahedronVolume(ndCoordList,centroid):
    """
    calculate the volume of a hexahedron from the coordinates of its vertices and its
    centroid
    
    The hexaherdon is subdivided into tetrahedons using the centroid of the
    hexahedron. Sub-tetrahedron volumes are calculated with
    :func:`tetrahedronVolume` and summed up.
    
    :param ndCoordList: list of arrays with node coordinates
    :type ndCoordList:  numpy.array
    :param centroid:    array with hexahedron centroid
    :type centroid:     numpy.array
    :return:            volume of tetrahedron
    :type:              float
    """
    vol = 0
    for nodeOrder in ([[0,1,2],[2,3,0],[4,5,6],[6,7,4],[0,4,7],[7,3,0],
    [1,5,6],[6,2,1],[2,3,7],[7,6,2],[1,0,4],[4,5,1]]):
        tetraCoordList = []
        tetraCoordList.append(ndCoordList[nodeOrder[0]])
        tetraCoordList.append(ndCoordList[nodeOrder[1]])
        tetraCoordList.append(ndCoordList[nodeOrder[2]])
        tetraCoordList.append(centroid)
        vol += tetrahedronVolume(tetraCoordList)
    return vol

def tetrahedron_volume_vert(vertList):
    """
    calculate the volume of a tetrahedron from its vertices
    
    :param vertList:          list of vertex objects
    :type vertList:           :class:`pyCFD_mesh.vertex.Vertex`
    :return:                  volume of tetrahedron
    :type:                    float
    """
    pointACoords = vertList[0].coords
    pointBCoords = vertList[1].coords
    pointCCoords = vertList[2].coords
    pointDCoords = vertList[3].coords
    vectDA = numpy.add(pointACoords,-pointDCoords)
    vectDB = numpy.add(pointBCoords,-pointDCoords)
    vectDC = numpy.add(pointCCoords,-pointDCoords)
    volume = math.fabs(numpy.dot(vectDA,numpy.cross(vectDB,vectDC)) / 6)
    return volume

def prism_volume_vert(vertList):
    """
    calculate the volume of a prism from its vertices
    
    The prism is subdivided into tetrahedons. Sub-tetrahedron volumes are
    calculated with :func:`tetrahedron_volume_vert` and summed up.
    
    :param vertList:          list of vertex objects
    :type vertList:           :class:`pyCFD_mesh.vertex.Vertex`
    :return:                  volume of prism
    :type:                    float
    """
    vol = 0
    ndCoordList = []
    for vertI in vertList:
        ndCoordList.append(vertI.get_coords)
    for nodeOrder in ([[0,1,2],[2,3,0],[2,5,3],[1,4,3],[1,3,0],[3,4,5]]):
        tetraCoordList = []
        tetraCoordList.append(ndCoordList[nodeOrder[0]])
        tetraCoordList.append(ndCoordList[nodeOrder[1]])
        tetraCoordList.append(ndCoordList[nodeOrder[2]])
        tetraCoordList.append(vertList[-1])
        vol += tetrahedron_volume_vert(tetraCoordList)
    return vol

def hexahedron_volume_vert(vertList):
    """
    calculate the volume of a hexahedron from its vertices
    
    The hexahedron is subdivided into tetrahedons. Sub-tetrahedron volumes are
    calculated with :func:`tetrahedron_volume_vert` and summed up.
    
    :param vertList:          list of vertex objects
    :type vertList:           :class:`pyCFD_mesh.vertex.Vertex`
    :return:                  volume of hexahedron
    :type:                    float
    """
    vol = 0
    for nodeOrder in ([[0,1,2],[2,3,0],[4,5,6],[6,7,4],[0,4,7],[7,3,0],
    [1,5,6],[6,2,1],[2,3,7],[7,6,2],[1,0,4],[4,5,1]]):
        tetraCoordList = []
        tetraCoordList.append(vertList[nodeOrder[0]])
        tetraCoordList.append(vertList[nodeOrder[1]])
        tetraCoordList.append(vertList[nodeOrder[2]])
        tetraCoordList.append(vertList[-1])
        vol += tetrahedron_volume_vert(tetraCoordList)
    return vol
   
def polyhedron_centroid_and_volume(cell, face_list):
    """
    calculate the centroid and volume of a polygon from the list of its faces
    
    **Calculation steps**
    
    * polyhedra is subdivided into pyramids using the geometric centroid of the polyhedra.
      Therefore the cell itself is needed as input with its vertices updated.

    * centers of the sub-pyramids and their volume is calculated

    * polyhedra centroid is calculated as volume weighted average of the sub pyramids

    * polyhedra volume is calculated as the sum of sub pyramids volumes
    
    :param cell:      incomplete cell object
    :type cell:       pyCFD_mesh.cell.Cell
    :param face_list: list of face objects defining the cell
    :type face_list:  dict
    :return:          dictionary [weighted centroid, volume]
    :rtype:           dict
    """
    number_of_faces = len(face_list)
    number_of_vertices = len(cell.vertices)
    if number_of_faces < 4:
        sys.exit("invalid cell element with low number of vertices (< 4)!")
    vertex_coords = numpy.array([cell.vertices[i_vertex].coords for i_vertex in xrange(number_of_vertices)])
    centroid = numpy.array([vertex_coords.sum(0) / float(number_of_vertices),] * number_of_faces)
    centroid_0 = centroid[0]
    face_centers = numpy.array([face_list[i_face].C for i_face in xrange(number_of_faces)])
    pyramid_centers = (0.75*face_centers + 0.25*centroid)
    vector_apex_to_base = numpy.add(face_centers, -centroid)
    Sf = numpy.array([face_area_vect_vert(face_list[i_face].vertices,centroid_0) for i_face in xrange(number_of_faces)])
    pyramid_volumes = numpy.einsum('ij,ij->i',vector_apex_to_base, Sf) / 3.
    # calculate poligon weighted centroid and area    
    weighted_centroid = numpy.array([pyramid_centers[i] * pyramid_volumes[i] for i in xrange(number_of_faces)])
    vol_sum = pyramid_volumes.sum()
    weighted_centroid_sum = weighted_centroid.sum(0) / vol_sum
    result = [weighted_centroid_sum, vol_sum]
    return result   






