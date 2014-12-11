"""
module for mesh objects (cells/faces)
"""
__author__ = "Bence Somogyi"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Bence Somogyi"
__email__ = "bencesomogyi@ivt.tugraz.at"
__status__ = "Prototype"

class MeshObject:
    """
    abstract class for mesh objects with father and vertices
    """
    def __init__(self):
        self.vertices = []
        """list of vertices defining the object"""
        self.father = []
        """reference to father object"""
        
def are_mesh_objects_equal(obj1, obj2):
    """
    compare mesh cells or faces for equality
    """
    import pyCFD_mesh.face as face
    import pyCFD_mesh.cell as cell
    import pyCFD_mesh.vertex as vertex
    # find out mode: cells/faces
    if isinstance(obj1, face.Face):
        face_mode = True
    if isinstance(obj1, cell.Cell):
        face_mode = False
    cell_mode = not face_mode
    
    # check if inputs are the same type of objects
    if face_mode:
        if (isinstance(obj1, face.Face) is not True):
            print 'use the same instances for this function!'
            return False
    if cell_mode:
        if (isinstance(obj1, cell.Cell) is not True):
            print 'use the same instances for this function!'
            return False
    vertex_match = 0
    
    # check if objs has the same # of vertices
    if (len(obj1.vertices) != len(obj2.vertices)):
        return False

    # check if obj vertices are matching
    check = True
    for vert1 in obj1.vertices:
        for vert2 in obj2.vertices:
            if vertex.are_vertices_equal(vert1, vert2):
                vertex_match += 1
    if vertex_match != len(obj1.vertices):
        check = False
    return check
