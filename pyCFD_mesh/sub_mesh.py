import pyCFD_mesh.generic_mesh as mesh

class SubMesh(mesh.GenericMesh):
    """
    Class for submeshes with a given list of cells and faces (e.g only the
    internal faces of a mesh)
    """
    def __init__(self,cellList=None,faceList=None):
        """
        constructor for the SubMesh class
        """
        mesh.GenericMesh.__init__(self)
        """
        **Constructor**
        
        :param cellList: default: None, list of vertices defining the face
        :type cellList:  :class:`pyCFD_mesh.cell.Cell`
        :param faceList: default: None, list of vertices defining the face
        :type faceList:  :class:`pyCFD_mesh.cell.Cell`
        """
        
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