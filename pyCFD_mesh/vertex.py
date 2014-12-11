#-------------------------------------------------------------------------------
# Name:        vertex
# Purpose:     class for 3D vertices
#
# Author:      bencesomogyi
#
# Created:     19.10.2013
# Copyright:   (c) bencesomogyi 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import numpy

class Vertex:
    """
    class for 3D vertices
    """
    def __init__(self,X=0.0,Y=0.0,Z=0.0):
        """
        **Constructor**
        
        :param X: default: 0.0, x coordinate of the vertex
        :type X:  float
        :param Y: default: 0.0, y coordinate of the vertex
        :type Y:  float
        :param Z: default: 0.0, z coordinate of the vertex
        :type Z:  float
        """
        self.X = X
        """vertex X coordinate"""
        self.Y = Y
        """vertex Y coordinate"""
        self.Z = Z
        """vertex Z coordinate"""
        self.coords = numpy.array([X, Y, Z])
        """vector of vertex coordinates"""
        self.father = []
        """reference to father object"""
        self.faces = []
        """reference to connected faces"""        
        self.cells = []
        """reference to connected cells"""
        self.id = 0.
        """vertex id"""

    def get_coords(self):
        """return coordinates as numpy array"""
        return numpy.array([self.X,self.Y,self.Z])

    def setX(self,newX):
        """setter for X coordinate of vertex"""
        self.X = newX
        self.coords[0] = newX

    def setY(self,newY):
        """setter for Y coordinate of vertex"""
        self.Y = newY
        self.coords[1] = newY

    def setZ(self,newZ):
        """setter for Z coordinate of vertex"""
        self.Z = newZ
        self.coords[2] = newZ
        
    def print_coordinates(self):
        print "id: "+str(self.id)+" "+str(self.X)+" "+str(self.Y)+" "+str(self.Z)
        
    def get_cell_ids(self):
        cell_ids = []
        for cell_ in self.cells:
            cell_ids.append(cell_.id)
        return cell_ids

#    def __hash__(self):
#        return self.X+10*self.Y+100*self.Z

    def __eq__(self,other):
        return (self.X == other.X) and (self.Y == other.Y) and (self.Z == other.Z)

def are_vertices_equal(vertex1,vertex2):
    return (vertex1.X == vertex2.X) and (vertex1.Y == vertex2.Y) and (vertex1.Z == vertex2.Z)
    
def get_independent_vertices(vertex_list):
    """
    return the list of independent vertices
    
    :param vertex_list: list of vertices
    :type vertex_list: dict
    :return: list of independent vertices
    :rtype: dict
    """
    indep_vertices = []
    for vertex in vertex_list:
        if len(indep_vertices) == 0:
            indep_vertices.append(vertex)
        match = False
        for indep_vertex in indep_vertices:
            if vertex == indep_vertex:
                match = True
                break
        if match:
            continue
        else:
            indep_vertices.append(vertex)
    return indep_vertices
    
def get_list_of_ids(vertex_list):
    """
    return the id / list of ids for a vertex / list of vertices
    """
    id_list = []
    if isinstance(vertex_list, Vertex):
        id_list.append(vertex_list.id)
    else:
        for vertex_ in vertex_list:
            id_list.append(vertex_.id)
    return id_list
            