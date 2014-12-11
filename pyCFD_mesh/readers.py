"""
module for reading existing meshes
"""
__author__ = "Bence Somogyi"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Bence Somogyi"
__email__ = "bencesomogyi@ivt.tugraz.at"
__status__ = "Prototype"

import sys
import time
import linecache
import pyCFD_mesh.vertex as vertex
import pyCFD_mesh.face as face
import pyCFD_mesh.cell as cell
import pyCFD_mesh.patch as patch
import pyCFD_mesh.generic_mesh as mesh

class MSHMesh(mesh.GenericMesh):
    """
    A class for reading gmsh meshes. The reader supports hex, tetra and wedge
    cells. Mesh data is read from _MESH/cells.msh and _MESH/boundary.msh.
    """
    def __init__(self):
        t_start = time.time()
        print ''
        print 'Processing mesh...'
        mesh.GenericMesh.__init__(self)
        
        self.FLAG_verbose = False
        """true/false to plot extra info when constructing
        the mesh > for debugging"""
        
        meshFile = '_MESH/cells.msh'
        f = open(meshFile, 'r')
        i = 1
        for line in f:
            if line.strip() == '$PhysicalNames':
                ptch_sl = i
            elif line.strip() == '$Nodes':
                nd_sl = i
            elif line.strip() == '$Elements':
                elm_sl = i
            i += 1
        f.close()

        # check for existence of patches, elements and nodes
        try:
            ptch_sl
        except NameError:
            print "no patch defined in msh file: ", meshFile
            # sys.exit("exit script...")
            ptchNo = 0
        else:
            ptchNo = int(linecache.getline(meshFile, ptch_sl+1))

        try:
            nd_sl
        except NameError:
            print "no nodes defined in cell.msh file!"
            sys.exit("exit script...")
        else:
            ndNo   = int(linecache.getline(meshFile, nd_sl+1))

        try:
            elm_sl
        except NameError:
            print "no elements defined in msh file!"
            sys.exit("exit script...")
        else:
            elmNo  = int(linecache.getline(meshFile, elm_sl+1))

        # fill data lists
        #   Loop patches and fill list with patch names
        patchNames = []
        for i in xrange(ptchNo):
            current_name = linecache.getline(meshFile, ptch_sl+2+i).split()[2]
            patchNames.append(current_name[1:-1])
        self.patchNames = patchNames
        """list of patch names"""
        self.patches = []
        """list of patch objects"""

        #   loop nodes and fill list of nodes
        vertexList = []
        for i in xrange(ndNo):
            new_node = vertex.Vertex()
            new_node.setX(float(linecache.getline(meshFile, nd_sl+2+i).split()[1]))
            new_node.setY(float(linecache.getline(meshFile, nd_sl+2+i).split()[2]))
            new_node.setZ(float(linecache.getline(meshFile, nd_sl+2+i).split()[3]))
            new_node.father.append(self)
            new_node.id = i
            vertexList.append(new_node)
        self.vertices = vertexList
        """list of vertex objects"""

        #   loop elements and fill lists for element type, node list per element
        if self.FLAG_verbose:        
            print 'creating cells...'
            t_0 = time.time()
        elementList = []
        self.faces = []
        """list of independent face objects"""
        for i in xrange(elmNo):
            actual_type = int(linecache.getline(meshFile, elm_sl+2+i).split()[1])
            actual_node_number = 0
            if ((actual_type != 4) and (actual_type != 5) and (actual_type != 6)):
                continue
            if actual_type == 4:
                # 4-node tetrahedron
                actual_node_number = 4
            elif actual_type == 5:
                # 8-node hexahedron
                actual_node_number = 8
            elif actual_type == 6:
                # 6-node prism
                actual_node_number = 6
            actual_vertices = []
            for j in xrange(actual_node_number):
                node_i = int(linecache.getline(meshFile, elm_sl+2+i).split()[5+j])-1
                actual_vertices.append(vertexList[node_i])
            temp_cell = cell.Cell(actual_vertices)
            temp_cell.father.append(self)
            elementList.append(temp_cell)
            temp_cell.id = i
        self.cells = elementList
        """list of cell objects"""
        self.faces = filter(None, self.faces)
        for face_i in self.faces:
            face_i.father.append(self)
        if self.FLAG_verbose:
            print 'number of faces: '+str(len(self.faces))
            print 'DONE in '+str(time.time()-t_0)+' s'

        #   find boundary faces, update Sf, weights and id
        if self.FLAG_verbose:
            print 'update "isBnd" flag and surface vectors...'
            t = time.time()
        for i, face_i in enumerate(self.faces):
            if (len(face_i.cells) == 1):
                face_i.isBnd = True
            # update face normal vectors
            face_i.update_Sf()
            # update face weights
            face_i.update_weights()
            # update face id
            face_i.id = i
        if self.FLAG_verbose:
            print 'DONE in '+str(time.time()-t)+' s'

        self.meshFile = meshFile
        """Name of mesh file loaded."""

        self.get_boundary_info()
        
        print '   DONE in '+str(time.time()-t_start)+' s'

    def get_boundary_info(self):
        """
        read boundary information from _MESH/bounday.msh file
        """
        #meshName,meshExt = os.path.splitext(self.meshFile)
        bnd_mesh_file = '_MESH/boundary.msh'
        f = open(bnd_mesh_file, 'r')
        i = 1
        for line in f:
            if line.strip() == '$Elements':
                elm_sl = i
            elif line.strip() == '$Nodes':
                nd_sl = i
            i += 1
        f.close()

        try:
            elm_sl
        except NameError:
            print "no elements defined in boundary.msh file!"
            sys.exit("exit script...")
        else:
            elmNo  = int(linecache.getline(bnd_mesh_file, elm_sl+1))
            #int(linecache.getline(bndMeshFile, elm_sl+1))
        try:
            nd_sl
        except NameError:
            print "no nodes defined in boundary.msh file!"
            sys.exit("exit script...")
        else:
            ndNo   = int(linecache.getline(bnd_mesh_file, nd_sl+1))

        #   loop nodes and fill list of nodes
        bnd_vertex_list = []
        if self.FLAG_verbose:        
            print 'reading boundary vertices...'
            t_0 = time.time()
        for i in xrange(ndNo):
            newNode = vertex.Vertex()
            newNode.setX(float(linecache.getline(bnd_mesh_file, nd_sl+2+i).split()[1]))
            newNode.setY(float(linecache.getline(bnd_mesh_file, nd_sl+2+i).split()[2]))
            newNode.setZ(float(linecache.getline(bnd_mesh_file, nd_sl+2+i).split()[3]))
            newNode.father.append(self)
            bnd_vertex_list.append(newNode)
        if self.FLAG_verbose:
            print 'DONE in '+str(time.time()-t_0)+' s'

        # exchange bnd vertices with existing vertices
        if self.FLAG_verbose:
            print 'exchange boundary vertices with cell vertices...'
            t_0 = time.time()
        for ind_bnd_vertex, bnd_vertex in enumerate(bnd_vertex_list):
            indMeshVert = self.vertices.index(bnd_vertex)
            bnd_vertex_list.pop(ind_bnd_vertex)
            bnd_vertex_list.insert(ind_bnd_vertex, self.vertices[indMeshVert])
        if self.FLAG_verbose:
            print 'DONE in '+str(time.time()-t_0)+' s'

        # update bndIds for the faces of the mesh
        if self.FLAG_verbose:
            print 'find bndIds for faces of mesh and add to patch...'
            t_0 = time.time()
        #  create empty patches
        for patchI in xrange(len(self.patchNames)):
            self.patches.append(patch.Patch([], self.patchNames[patchI]))
        #  find bndIds
        for faceI in xrange(elmNo):
            patchI = int(linecache.getline(bnd_mesh_file, elm_sl+2+faceI).split()[3])-1
            patchFaceVertices = []
            elementType = int(linecache.getline(bnd_mesh_file, elm_sl+2+faceI).split()[1])
            nodes = 0
            if elementType == 2:
                # 3-node trirangle
                nodes = 3
            elif elementType == 3:
                # 4-node quadrangle
                nodes = 4
            for j in xrange(nodes):
                patchFaceVertices.append(bnd_vertex_list[int(linecache.getline(bnd_mesh_file,elm_sl+2+faceI).split()[5+j])-1])
            tempFace = face.Face(patchFaceVertices)
            for faceI in patchFaceVertices[0].faces:
                if faceI == tempFace:
                    faceI.bndId = patchI
                    faceI.isBnd = True
                    self.patches[patchI].faces.append(faceI)
                    break
        if self.FLAG_verbose:
            print 'DONE in '+str(time.time()-t_0)+' s'
            
class FoamMesh(mesh.GenericMesh):
    """
    A class for reading OpenFOAM meshes. The reader supports hex, tetra and
    wedge cells. Mesh data is read from _MESH/dir/cells.msh and
    _MESH/dir/boundary.msh.
    """
    def __init__(self, dir=""):
        """
        **Constructor**
        
        :param dir: subdir within the _MESH directory
        :type dir:  string
        """
        t_start = time.time()
        print ''
        print 'Processing mesh...'
        mesh.GenericMesh.__init__(self)
        
        self.FLAG_verbose = False
        """true/false to plot extra info when constructing
        the mesh > for debugging"""
        
        if dir == "":
            path_ = "_MESH/"
        else:
            path_ = "_MESH/"+dir+"/"
        
        # ========        
        # VERTICES
        # ========
        
        
        if self.FLAG_verbose:
            print 'creating vertex objects...'
            t = time.time()
        
        # find nodes
        node_file = path_+'points'
        f = open(node_file, 'r')
        i = 1
        for line in f:
            if line.strip() == '(':
                nd_sl = i+1
                break
            i += 1
        f.close()
        try:
            nd_sl
        except NameError:
            print "node definition not found, check for formatting:"
            print "..."
            print "// * * * *"
            print ""
            print "# of nodes"
            print "( >>> THIS SHOULD BE A SEPARATE LINE!"
            print "(coordX coordY coordZ)"
            print "..."
            sys.exit("exit script...")
        else:
            ndNo   = int(linecache.getline(node_file, nd_sl-2))
            
        # create nodes
        vertexList = []
        for i_node in xrange(ndNo):
            new_node = vertex.Vertex()
            # get coords as a string list
            string = linecache.getline(node_file, nd_sl+i_node).split()
            string_new = []
            for str_elem in string:
                str_elem_list = list(str_elem)
                for i,char in enumerate(str_elem_list):
                    if char == "(" or char == ")":
                        str_elem_list[i]=""
                string_new.append("".join(str_elem_list))
            new_node.setX(float(string_new[0]))
            new_node.setY(float(string_new[1]))
            new_node.setZ(float(string_new[2]))
            new_node.father.append(self)
            new_node.id = i_node
            vertexList.append(new_node)
        self.vertices = vertexList
        """list of vertex objects"""
        if self.FLAG_verbose:
            print 'DONE in '+str(time.time()-t)+' s'
        
        # =====        
        # FACES
        # =====
        
        if self.FLAG_verbose:
            print 'creating face objects...'
            t = time.time()
        
        # find faces
        face_file = path_+'faces'
        f = open(node_file, 'r')
        i = 1
        for line in f:
            if line.strip() == '(':
                fc_sl = i+1
                break
            i += 1
        f.close()
        try:
            fc_sl
        except NameError:
            print "face definition not found, check for formatting:"
            print "..."
            print "// * * * *"
            print ""
            print "# of faces"
            print "( >>> THIS SHOULD BE A SEPARATE LINE!"
            print "# of vertices(vertex0 vertex2 ...)"
            print "..."
            sys.exit("exit script...")
        else:
            fcNo   = int(linecache.getline(face_file, fc_sl-2))
        
        # create faces
        face_list = []
        for i_face in range(fcNo):
            # get indices as a string list
            string = linecache.getline(face_file, fc_sl+i_face).split()
            string_new = []
            for str_elem in string:
                str_elem_list = list(str_elem)
                for i,char in enumerate(str_elem_list):
                    if len(str_elem_list) == 1:
                        continue
                    if (i == 0 and str_elem_list[1] == "(")  or char == "(" or char == ")":
                        str_elem_list[i]=""
                string_new.append("".join(str_elem_list))
            vertex_list = []
            for i_vertex in string_new:
                vertex_list.append(self.vertices[int(i_vertex)])
            new_face = face.Face(vertex_list)
            new_face.father.append(self)
            new_face.id = i_face
            face_list.append(new_face)
            for vertex_ in new_face.vertices:
                vertex_.faces.append(new_face)
        self.faces = face_list
        """list of face objects"""
        if self.FLAG_verbose:
            print 'DONE in '+str(time.time()-t)+' s'

        # =====        
        # CELLS
        # =====  

        if self.FLAG_verbose:
            print 'creating cell objects...'
            t = time.time()      
        
        # find cells
        owners_file = path_+'owner'
        neighbours_file = path_+'neighbour'
        f = open(owners_file, 'r')
        i = 1
        for line in f:
            if line.strip() == '(':
                o_sl = i+1
                break
            i += 1
        f.close()
        try:
            o_sl
        except NameError:
            print "owner definition not found, check for formatting:"
            print "..."
            print "// * * * *"
            print ""
            print "# of owners"
            print "( >>> THIS SHOULD BE A SEPARATE LINE!"
            print "face id"
            print "..."
            sys.exit("exit script...")
        else:
            face_number   = int(linecache.getline(owners_file, o_sl-2))
            if face_number != len(self.faces):
                sys.exit("inconsistent number of faces in _MESH/faces and _MESH/owner!")
        f = open(neighbours_file, 'r')
        i = 1
        for line in f:
            if line.strip() == '(':
                n_sl = i+1
                break
            i += 1
        f.close()
        try:
            n_sl
        except NameError:
            print "neighbour definition not found, check for formatting:"
            print "..."
            print "// * * * *"
            print ""
            print "# of neighbours"
            print "( >>> THIS SHOULD BE A SEPARATE LINE!"
            print "face id"
            print "..."
            sys.exit("exit script...")
        else:
            internal_face_number  = int(linecache.getline(neighbours_file, n_sl-2))
            
        # create cells
        owner_list = []
        neighbour_list = []
        internal_face_number = int(linecache.getline(neighbours_file, n_sl-2))
        faces_of_cells = []
        for i_face in range(face_number):
            owner_list.append(int(linecache.getline(owners_file, o_sl+i_face)))
        for i_face in range(internal_face_number):
            neighbour_list.append(int(linecache.getline(neighbours_file, n_sl+i_face)))
        cell_number = max(owner_list)+1
        for i_cell in range(cell_number):
            faces_of_cells.append([])
        for i_face in range(face_number):
            faces_of_cells[int(linecache.getline(owners_file, o_sl+i_face))].append(i_face)
        for i_face in range(internal_face_number):
            faces_of_cells[int(linecache.getline(neighbours_file, n_sl+i_face))].append(i_face)
        #  append faces to the cell's list
        for i_cell in range(cell_number):
            face_list = []
            for i_face in faces_of_cells[i_cell]:
                face_list.append(self.faces[i_face])
            temp_cell = cell.Cell(face_list)
            temp_cell.id = i_cell
            temp_cell.father.append(self)
            self.cells.append(temp_cell)
        del i_face, i_cell, o_sl, n_sl, face_list, faces_of_cells, owners_file, neighbours_file
        
        # append cells to vertices and faces
        for cell_ in self.cells:
            for vertex_ in cell_.vertices:
                vertex_.cells.append(cell_)
                
        # loop owners list
        for i_face in range(face_number):
            self.faces[i_face].cells.append(self.cells[owner_list[i_face]])
        
        # loop neighbours list
        for i_face in range(internal_face_number):
            self.faces[i_face].cells.append(self.cells[neighbour_list[i_face]])

        # update Sf, weights and isBnd
        for face_ in self.faces:
            face_.update_Sf()
            face_.update_weights()
            face_.update_gradient_weights()
            if len(face_.cells) == 1:
                face_.isBnd = True
                
        if self.FLAG_verbose:
            print 'DONE in '+str(time.time()-t)+' s'
                
        # =======        
        # PATCHES
        # =======
                
        if self.FLAG_verbose:
            print 'creating patches...'
            t = time.time()
        
        # find patches
        patches_file = path_+'boundary'
        f = open(patches_file, 'r')
        i = 1
        p_start = []
        p_stop  = []
        p_num_found = False
        for line in f:
            if line.strip() == '(':
                p_num = int(linecache.getline(patches_file, i-1))
                p_num_found = True
            if p_num_found and line.strip() == '{':
                p_start.append(i+1)
            if p_num_found and line.strip() == '}':
                p_stop.append(i)
            i += 1
        f.close()
        if len(p_stop) is not p_num:
            print "p_stop: "+str(len(p_stop))
            sys.exit("could not find all patches, stopping")
        elif len(p_start) is not p_num:
            print "p_start: "+str(len(p_start))
            sys.exit("could not find all patches, stopping")
        self.patches = [None] * p_num
        for i_patch in range(p_num):
            patch_name = str(linecache.getline(patches_file, p_start[i_patch]-2)).strip()
            for line_ in range(p_start[i_patch], p_stop[i_patch]):
                if linecache.getline(patches_file, line_).strip()[0:6] == "nFaces":
                    n_faces = int(linecache.getline(patches_file, line_).strip()[6:-1])
                if linecache.getline(patches_file, line_).strip()[0:9] == "startFace":
                    start_face = int(linecache.getline(patches_file, line_).strip()[9:-1])
            patch_face_list = []
            for i_face,face_ in enumerate(self.faces[start_face:start_face+n_faces]):
                patch_face_list.append(face_)
                face_.bndId = i_patch
                face_.inPatchId = i_face
            self.patches[i_patch] = patch.Patch(patch_face_list, patch_name)
            self.patchNames.append(patch_name)
        
        if self.FLAG_verbose:
            print 'DONE in '+str(time.time()-t)+' s'
            
        print 'MESH READ in '+str(time.time()-t_start)+' s'