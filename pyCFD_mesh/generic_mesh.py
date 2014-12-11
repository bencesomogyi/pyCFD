"""
module for generic meshes
"""
__author__ = "Bence Somogyi"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Bence Somogyi"
__email__ = "bencesomogyi@ivt.tugraz.at"
__status__ = "Prototype"

import sys
import numpy

class GenericMesh:
    """
    basic class for meshes
    """
    
    def __init__(self):
        # attributes
        self.vertices = []
        """List of all vertices"""
        self.faces = []
        """List of independent faces"""
        self.cells = []
        """List of cells in the mesh"""
        self.patchNames = []
        """List of patch names"""
        self.patches = []
        """List of patches"""
        self.fields = []
        """List of fields"""
        
    def plot_mesh_data(self):
        """
        plots basic mesh data
        
        * number of vertices
        
        * number of faces
        
            * number of boundary faces
            
            * number of internal faces
            
        * number of cells
        
        * number of patches
        
            * list of patches
        """
        print ''
        print '========================='
        print 'MESH DATA:'
        print '- number of vertices: '+str(len(self.vertices))
        print '- number of faces: '+str(len(self.faces))
        bndFaces = 0
        intFaces = 0
        for face in self.faces:
            if face != None:
                if face.isBnd == True:
                    bndFaces += 1
                else:
                    intFaces += 1
        print '-- number of boundary faces: '+str(bndFaces)
        print '-- number of internal faces: '+str(intFaces)
        print '- number of cells: '+str(len(self.cells))
        print '- number of patches: '+str(len(self.patches))
        for patchName in self.patchNames:
            print '-- '+str(patchName)
        print '========================='
        
    def get_patch(self, patch_name):
        """
        return the patch with matching name
        
        :param patch_name: name of patch
        :type patch_name:  string
        :return:           geometric patch object
        :rtype:            :class:`pyCFD_mesh.patch.Patch`
        """
        patch_found = -1
        for patch_i,patch_name_i in enumerate(self.patchNames):
            if patch_name == patch_name_i:
                patch_found = patch_i
                break
        
        if patch_found == -1:
            print patch_name + " does not exist in the mesh, possible patch names are:"
            for patch_name_i in self.patchNames:
                print patch_name_i
            sys.exit()
        
        return self.patches[patch_found]
        
    def get_volumes(self):
        """
        return volumes of all cells within the mesh as vector
        
        :return: array with cell volumes
        :rtype:  numpy.array
        """
        cell_volumes = numpy.zeros((len(self.cells),1))
        for i_ in range(len(cell_volumes)):
            cell_volumes[i_] = self.cells[i_].V
            
        return cell_volumes
        
    def get_areas(self):
        """
        return areas of all faces within the mesh as vector
        
        :return: array with face areas
        :rtype:  numpy.array
        """
        face_areas = numpy.zeros((len(self.faces),1))
        for i_ in range(len(face_areas)):
            face_areas[i_] = self.faces[i_].A
            
        return face_areas