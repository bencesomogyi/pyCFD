"""
module for cell faces
"""
__author__ = "Bence Somogyi"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Bence Somogyi"
__email__ = "bencesomogyi@ivt.tugraz.at"
__status__ = "Prototype"

import numpy
import sys
import pyCFD_mesh.vertex as vertex
import pyCFD_mesh.mesh_object as mesh_object
import pyCFD_geometric_tools.geomTools as geomTools

if sys.platform == 'win32':
    import pyCFD_operators.cython_boost_win32.cy_operators as cy_operators
elif sys.platform == 'linux2':
    import pyCFD_operators.cython_boost_linux2.cy_operators as cy_operators
else:
    sys.exit("unknown platform " + sys.platform + " found in geomTools.py, stopping...")

class Face(mesh_object.MeshObject):
    """class for cell faces"""
    def __init__(self,vertexList):
        """
        **Constructor**
        
        :param vertexList: list of vertices defining the face
        :type vertexList:  :class:`pyCFD_mesh.vertex.Vertex`
        """
        mesh_object.MeshObject.__init__(self)
        # check vertex input
        node_number = len(vertexList)
        if ((node_number != 3) and (node_number != 4)):
            print self, ' created with an invalid number of vertices, should be 3 or 4'
            return
        nodesIndependent = True
        for indVertI in range(len(vertexList)-1):
            for vertJ in vertexList[indVertI+1:]:
                if vertex.are_vertices_equal(vertexList[indVertI],vertJ):
                    nodesIndependent = False
        if nodesIndependent == False:
            print 'nodes are not independent in ', self
            sys.exit()
            return
            
        # get face centroid and area
        centr_and_area = geomTools.polygon_centroid_and_area(vertexList)
        
        # initialize Sf with 0,0,0
        Sf = numpy.array([0., 0., 0.])

        self.A = centr_and_area[1]
        """face areas [m2]"""
        self.Sf = Sf
        """face area vectors [m2]"""
        self.vertices = vertexList
        """list of vertices defining the face"""
        self.cells = []
        """list of connected cells"""
        self.isBnd = False
        """bool for boundary faces"""
        self.bndId = -1
        """boundary ID of face"""
        self.C = centr_and_area[0]
        """face centroid"""
        self.weights = []
        """weights for linear interpolation [w_owner, w_neigbour]"""
        self.gradWeights = []
        """weights for Gauss gradient correction [w_owner, w_neigbour]"""        
        self.ffToF = numpy.array([0., 0., 0.])
        """is the shortest vector starting from the line between owner and neighbour pointing to the face centroid"""
        self.id = 0
        """face id"""
        self.inPatchId = 0
        """face id in patch face list"""

    def __eq__(self,otherFace):
        return are_faces_equal(self,otherFace)
        
    def update_Sf(self):
        """
        update surface vectors
        """
        node_number = len(self.vertices)
        if ((node_number != 3) and (node_number != 4)):
            print self, ' created with an invalid number of vertices, should be 3 or 4'
            return
        # get face area
        Sf = 0
        if node_number == 3: # triangle
            Sf = geomTools.triangle_area_vect_vert(self.vertices,self.cells[0].C)
        elif node_number == 4: # quadrangle
            Sf = geomTools.quadrangle_area_vect_vert(self.vertices,self.cells[0].C)
        else:
            error_message = "update_Sf requested for face "+str(self.id)+" with "+str(node_number)+"vertices. Stopping..."
            sys.exit(error_message)
        self.Sf = Sf
        
    def update_weights(self):
        r"""
        update linear interpolation weights
        
        Calculation of the weights is done according to their normal distance
        from the face:
        
        .. image:: _images/update_linear_weights.png
            :width: 300px
            :align: center
        
        * | for the neighbour cell
          | :math:`g_N = \frac{\vec{d_{Of}} \cdot \vec{e_f}}{\vec{d_{Of}} \cdot \vec{e_f} + \vec{d_{fC}} \cdot \vec{e_f}}`
        
        * | for the owner cell
          | :math:`g_O = 1-g_N`
        
        , where :math:`\vec{e_f}` is the face normal unit vector
        :math:`\vec{e_f} = \frac{\vec{S_f}}{\|\vec{S_f}\|}` 
        """
        cell_num = len(self.cells)
        if ((cell_num != 2) and (cell_num != 1)):
            error_message = "wrong number of cells (" + str(cell_num) + ") for face" + str(self.id)
            sys.exit(error_message)
            
        if cell_num == 2:
            # cells
            #   owner cell
            cell_o = self.cells[0]
            #   neighbour cell
            cell_n = self.cells[1]
            
            # distances from plane
            #   d_Of dot e_f        
            of = numpy.dot(numpy.add(self.C, -cell_o.C), self.Sf/self.A)
            #   d_fN dot e_f
            fn = numpy.dot(numpy.add(cell_n.C, -self.C), self.Sf/self.A)
            
            # weights
            #   neighbour
            w_n = of / ( of + fn )
            #   owner
            w_o = 1.0 - w_n
            
            self.weights.append(w_o)
            self.weights.append(w_n)
            
    def update_gradient_weights(self):
        r"""
        update non-conjunctional interpolation location to be used in Gauss
        gradient iterations
        
        The vector :math:`\vec{r_{f'}}` belongs to the face and is calculated by:
            
        .. math::
            
            \vec{r_{f'}} = \vec{r_O} + \frac{\vec{r_{Of}} \cdot \vec{r_{ON}}}{\vec{r_{ON}} \cdot \vec{r_{ON}}} \left(\vec{r_O} - \vec{r_N} \right)
            
        .. image:: _images/update_grad_weights.png
            :width: 400px
            :align: center
        
        , :math:`\vec{r_{ff'}}` is chosen to be perpendicular to 
        :math:`\vec{r_{ON}}`.
        
        Interpolation weights based on this fictious point :math:`f'` are calculated the following way:
        
        * | for the owner:
          | :math:`g_O = \frac{\left| \vec{r_N} - \vec{r_{f'}} \right|}{\left| \vec{r_N} - \vec{r_O} \right|}`
          
        * | for the neighbour:
          | :math:`g_N = 1- g_O`
        """
        cell_num = len(self.cells)
        if ((cell_num != 2) and (cell_num != 1)):
            error_message = "wrong number of cells (" + str(cell_num) + ") for face" + str(self.id)
            sys.exit(error_message)
        
        if cell_num == 2:
            # cells
            #   owner cell
            cell_o = self.cells[0]
            #   neighbour cell
            cell_n = self.cells[1]
            
            ff = cell_o.C + numpy.dot(
                                      numpy.add(self.C, -cell_o.C),
                                      numpy.add(cell_n.C, -cell_o.C)
                                     ) \
                          / numpy.dot(
                                      numpy.add(cell_n.C, -cell_o.C),
                                      numpy.add(cell_n.C, -cell_o.C)
                                     ) \
                          * numpy.add(cell_n.C, -cell_o.C)
                
            # weights
            #   owner
            w_o = cy_operators.cy_linalg_norm(numpy.add(cell_n.C, -ff)) \
                / cy_operators.cy_linalg_norm(numpy.add(cell_n.C, -cell_o.C))
            #   neighbour
            w_n = 1.0 - w_o
            
            self.gradWeights.append(w_o)
            self.gradWeights.append(w_n)
            self.ffToF = numpy.add(self.C, -ff)
            
    def print_vertex_ids(self):
        """
        print the ids of the face's vertices.
        
        .. note:
            Vertex id is linked to a mesh.
        """
        print "face id: "+str(self.id)
        string = ""
        for vertex_ in self.vertices:
            string += str(vertex_.id)+" "
        print string
        
    def get_vertex_ids(self):
        """
        return the ids of the face's vertices.
        
        :return: list of vertex indices
        :rtype:  int
        
        .. note:
            Vertex id is linked to a mesh.
        """
        vertex_ids = []
        for vertex_ in self.vertices:
            vertex_ids.append(vertex_.id)
        return vertex_ids
        
    def get_cell_ids(self):
        """
        return the id of the cell which owns the face
        
        :return: id of the owner cell
        :rtype:  int
        """
        cell_ids = []
        for cell_ in self.cells:
            cell_ids.append(cell_.id)
        return cell_ids
        
    def get_Sf(self, cell_):
        """
        return face area vector pointing out from one of the cells that own the
        face
        
        :param cell_: reference cell
        :type cell_:  :class:`pyCFD_mesh.cell.Cell`
        :return:      array of face area vector
        :rtype:       numpy.array
        """
        if cell_.id == self.cells[0].id:
            return self.Sf
        elif cell_.id == self.cells[1].id:
            return self.Sf * -1.0
        else:
            sys.exit("wrong face.cells list found in Face.get_Sf, stopping!")
            
    def get_Sf_sign(self, cell_):
        """
        return
        
        * 1.0 if cell is the owner
        
        * -1.0 if cell is the neighbour
        
        :param cell_: reference cell
        :type cell_:  :class:`pyCFD_mesh.cell.Cell`
        :return:      1.0 or -1.0
        :rtype:       float
        """
        if cell_.id == self.cells[0].id:
            return  1.0
        elif cell_.id == self.cells[1].id:
            return -1.0
        else:
            sys.exit("wrong face.cells list found in Face.get_Sf, stopping!")
        
    def find_neighbour_vertex_on_face(self, vertex_):
        """
        returns a list of vertex objects which are the neighbours of reference
        vertex on the same face
        
        :param vertex_: reference vertex
        :type vertex_:  :class:`pyCFD_mesh.vertex.Vertex`
        :return:        list vertex objects
        :rtype:         :class:`pyCFD_mesh.vertex.Vertex`
        """
        neighbours = []
        n_ids = []
        if vertex_ not in self.vertices:
            sys.exit("vertex is not on the vertex list of face defined. Stopped in face.find_neighbour_vertex_on_face!")
        in_face_id = self.vertices.index(vertex_)
        if in_face_id == 0:
            n_ids.append(-1)
            n_ids.append(1)
        elif in_face_id == len(self.vertices)-1:
            n_ids.append(0)
            n_ids.append(-2)
        else:
            n_ids.append(in_face_id-1)
            n_ids.append(in_face_id+1)
        for v_id in n_ids:
            neighbours.append(self.vertices[v_id])
        return neighbours
        
def are_faces_equal(face1,face2):
    """
    check if two face objects are the same
    
    :param face1: first face to compare
    :type face1:  :class:`pyCFD_mesh.face.Face`
    :param face2: second face to compare
    :type face2:  :class:`pyCFD_mesh.face.Face`
    :return:      True or False
    :rtype:       bool
    """
    # check for None
    if face1 is None or face2 is None:
        return False

    # check if faces has the same # of vertices
    if ( len(face1.vertices) != len(face2.vertices) ):
        return False

    # check if face centroid is the same
    if ( not(numpy.array_equal(face1.C,face2.C)) ):
        return False
        
    # check if face normal is the same
    if ( not(numpy.allclose(face1.Sf,face2.Sf)) and not(numpy.allclose(face1.Sf,(-1.0*face2.Sf))) ):
        return False
        
    return True
