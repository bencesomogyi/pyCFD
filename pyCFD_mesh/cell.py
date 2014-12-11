"""
module for mesh cells
"""
__author__ = "Bence Somogyi"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Bence Somogyi"
__email__ = "bencesomogyi@ivt.tugraz.at"
__status__ = "Prototype"

import sys
import numpy
import pyCFD_mesh.vertex as vertex
import pyCFD_mesh.face as face
import pyCFD_mesh.mesh_object as mesh_object
import pyCFD_geometric_tools.geomTools as geomTools

class Cell(mesh_object.MeshObject):
    """class for cells"""
    def __init__(self,obj_list):
        """
        **Constructor**
        
        :param obj_list: list of objects defining the cell
        :type obj_list:  :class:`pyCFD_mesh.vertex.Vertex` or :class:`pyCFD_mesh.face.Face`
        """
        mesh_object.MeshObject.__init__(self)
        
        self.V = 0.0
        """cell volume [m3]"""
        self.vertices = []
        """list of vertices defining the cell"""
        self.faces = []
        """list faces of the cell"""
        self.C = numpy.array([0.0,0.0,0.0])
        """cell centroid coordinates"""
        self.id = 0
        """cell id"""
        
        if isinstance(obj_list[0], vertex.Vertex):
            self.create_cell_form_vertices(obj_list)
        elif isinstance(obj_list[0], face.Face):
            self.create_cell_from_faces(obj_list)
        else:
            raise TypeError(
            "unsupported input type '{}' in consturctor of '{}'"
            ).format(type(obj_list[0]), self.__class__)
        
    def create_cell_from_vertices(self, vertexList):
        """
        create cell from vertices
        
        :param vertexList: list of vertices
        :type vertexList:  :class:`pyCFD_mesh.vertex.Vertex`
        """
        # check vertex input
        nodeNumber = len(vertexList)
        if ((nodeNumber != 4) and (nodeNumber != 6) and (nodeNumber != 8)):
            print self, ' created with an invalid number of vertices, should be 4, 6 or 8'
            return
        nodesIndependent = True
        for indVertI in range(len(vertexList)-1):
            for vertJ in vertexList[indVertI+1:]:
                if vertex.are_vertices_equal(vertexList[indVertI],vertJ):
                    nodesIndependent = False
        if nodesIndependent == False:
            print 'nodes are not independent in  ', self
            return

        # get face list
        self.faces = []
        centroid = vertex.Vertex()
        V = 0.0
        if   nodeNumber == 4: # tetrahedron
            tetraCell = True
            prismCell = False
            hexaCell  = False
            fullOrder = [[0,1,2],[0,3,2],[3,1,2],[0,1,3]]
        elif nodeNumber == 6: # prism:
            tetraCell = False
            prismCell = True
            hexaCell  = False
            fullOrder = [[0,1,2],[3,4,5],[0,3,5,2],[0,1,4,3],[1,2,5,4]]
        elif nodeNumber == 8: # hexahedron:
            tetraCell = False
            prismCell = False
            hexaCell  = True
            fullOrder = [[0,1,2,3],[4,5,6,7],[0,1,5,4],[3,2,6,7],[0,3,7,4],[1,2,6,5]]
        for nodeOrder in fullOrder:
            vertexListCurr = []
            for i in nodeOrder:
                vertexListCurr.append(vertexList[i])
            tempFace = face.Face(vertexListCurr)
            for vertI in vertexListCurr:
                match = False
                for faceI in vertI.faces:
                    if faceI == tempFace:
                        match = True
                        break
                if match:
                    break
            if match: # existing face
                self.faces.append(faceI)
            else:     # new face
                self.faces.append(tempFace)
                for vertI in vertexListCurr:
                    vertI.faces.append(tempFace)
                vertexList[0].father[0].faces.append(tempFace)
        
        # get volume of cell
        if tetraCell:
            centroidCoords = geomTools.tetrahedron_centroid_vert(vertexList)
        elif prismCell:
            centroidCoords = geomTools.prism_centroid_vert(vertexList)
        elif hexaCell:
            centroidCoords = geomTools.hexahedron_centroid_vert(vertexList)
        centroid.setX(centroidCoords[0])
        centroid.setY(centroidCoords[1])
        centroid.setZ(centroidCoords[2])
        verticesWCentroid = []
        for vertexI in vertexList:
            verticesWCentroid.append(vertexI)
        verticesWCentroid.append(vertex.Vertex(centroid.X,centroid.Y,centroid.Z))
        if tetraCell:
            V = geomTools.tetrahedron_volume_vert(verticesWCentroid)
        elif prismCell:
            V = geomTools.prism_volume_vert(verticesWCentroid)
        elif hexaCell:
            V = geomTools.hexahedron_volume_vert(verticesWCentroid)

        # append cell to face list
        for faceI in self.faces:
            faceI.cells.append(self)

        self.V = V
        self.vertices = vertexList
        self.faces = self.faces
        self.C = centroidCoords
  
    def create_cell_from_faces(self, face_list):
        """
        create cell from faces
        
        :param vertexList: list of faces
        :type vertexList:  :class:`pyCFD_mesh.face.Face`
        """
        # find list of independent vertices
        vertex_list = []
        for face_ in face_list:
            for vertex_ in face_.vertices:
                vertex_list.append(vertex_)
        indep_vertices = vertex.get_independent_vertices(vertex_list)
        self.vertices.extend(indep_vertices)
        
        # add faces
        self.faces.extend(face_list)
        
        # calculate volume and centroid
        centr_and_vol = geomTools.polyhedron_centroid_and_volume(self,face_list)
        
        self.C = centr_and_vol[0]
        self.V = centr_and_vol[1]
        
    def print_face_ids(self):
        """
        print the ids of the cell's faces.
        
        .. note:
            Face id is linked to a mesh.
        """
        print "cell id: "+str(self.id)
        string = "faces: "
        for face_ in self.faces:
            string += str(face_.id)+" "
        print string
        
    def print_vertex_ids(self):
        """
        print the ids of the cell's vertices.
        
        .. note:
            Vertex id is linked to a mesh.
        """
        print "cell id: "+str(self.id)
        string = "vertices: "
        for vertex_ in self.vertices:
            string += str(vertex_.id)+" "
        print string
            

def set_father_of_cell_faces(cell):
    """
    set father for all faces of a cell
    
    :param cell: a cell object
    :type cell:  :class:`Cell`
    """
    if len(cell.father) != 1:
        sys.exit("cell must have 1 father!")

    for face in cell.faces:
        if len(face.father) == 0:
            face.father.append(cell.father[0])
        else:
            face.father[0] = cell.father[0]