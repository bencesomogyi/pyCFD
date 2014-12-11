#-------------------------------------------------------------------------------
# Name:        main
# Purpose:
#
# Author:      bencesomogyi
#
# Created:     20.02.2013
# Copyright:   (c) bencesomogyi 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------

#import os
#import numpy
#import time

#import pyCFD_mesh.vertex as vertex
#import pyCFD_mesh.cell as cell
#import pyCFD_mesh.face as face
import pyCFD_mesh.readers as readers
#import geomTools
#import pyCFD_VTK_tools.vtkTools as vtkTools
#import pyCFD_test_functions.testFunctions as testFunctions
import pyCFD_output.output as output
#import pyCFD_fields.fields as fields
#import pyCFD_fields.initialization as init

# create mesh object and read mesh data from MESHES directory
#myMesh = readers.MSHMesh()
myMesh = readers.FoamMesh("inclinedBlock")
myMesh.plot_mesh_data()
# create scalar field p and vector field U
#UInit = numpy.array([1.0,1.0,-1.0])
#p = fields.ScalarField(myMesh,'p',100)
#init.init_linear_scalar_sphere_distribution(p,100,-15,numpy.array([15.0,15.0,0.5]),30)
#U = fields.VectorField(myMesh,'U',UInit)
#output.write_mesh_file_with_fields(numpy.array([p,U]))
output.write_patch_files(myMesh)
#Sf = fields.SurfaceVectorField(myMesh, "Sf")
#Sf_vec = [myMesh.faces[i].Sf for i in xrange(len(myMesh.faces))]
#Sf.initialize_with_vector(Sf_vec)
#output.write_mesh_faces_with_field(Sf)
#output.write_mesh_file(myMesh)
#output.write_mesh_internal_faces(myMesh)
#output.write_mesh_faces(myMesh)

print ""
print "tetra not tested yet!"

print ""
print "testMesh FINSHED"