#-------------------------------------------------------------------------------
# Name:        genericMesh
# Purpose:     class for meshes defined by existing cells or faces
#
# Author:      bencesomogyi
#
# Created:     21.11.2013
# Copyright:   (c) bencesomogyi 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import pyCFD_mesh.generic_mesh as mesh

class Mesh(mesh.GenericMesh):
    """
    class for meshes defined by existing cells or faces
    """
    def __init__(self,cellList=None,faceList=None):
        """
        **Constructor**
        
        :param cellList: default: None, list of vertices defining the face
        :type cellList:  :class:`pyCFD_mesh.cell.Cell`
        :param faceList: default: None, list of vertices defining the face
        :type faceList:  :class:`pyCFD_mesh.cell.Cell`
        """
        mesh.GenericMesh.__init__(self)
        
        if cellList != None:
            if len(cellList) == 1:
                self.cells.append(cellList)
            else:
                self.cells.extend(cellList)
                
        if faceList != None:
            if len(faceList) == 1:
                self.faces.append(faceList)
            else:
                self.faces.extend(faceList)