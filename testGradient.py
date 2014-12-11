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

import numpy
#import math

import pyCFD_mesh.readers as readers
import pyCFD_output.output as output
import pyCFD_fields.fields as fields
import pyCFD_fields.calculated_fields as calculated_fields
#import pyCFD_operators.explicit_operators as exOp

# create mesh object and read mesh data from MESHES directory
myMesh = readers.FoamMesh()

# create and initialize field phi
phi = fields.ScalarField(myMesh, "phi")
for i_cell in range(len(phi.V)):
    cell_x = myMesh.cells[i_cell].C[0]
    cell_y = myMesh.cells[i_cell].C[1]
    phi.V[i_cell] =  pow(cell_x, 2.) * pow(cell_y, 2.)
for i_face in range(len(phi.A)):
    cell_x = myMesh.faces[i_face].C[0]
    cell_y = myMesh.faces[i_face].C[1]
    phi.A[i_face] =  pow(cell_x, 2.) * pow(cell_y, 2.)
    
# create and initialize field grad_phi_an
grad_phi_an = fields.VectorField(myMesh, "grad_phi_an")
for i_cell in range(len(phi.V)):
    cell_x = myMesh.cells[i_cell].C[0]
    cell_y = myMesh.cells[i_cell].C[1]
    grad_phi_an.V[i_cell][0] =  2*cell_x * pow(cell_y, 2.)
    grad_phi_an.V[i_cell][1] =  2*cell_y * pow(cell_x, 2.)
for i_face in range(len(phi.A)):
    cell_x = myMesh.faces[i_face].C[0]
    cell_y = myMesh.faces[i_face].C[1]
    grad_phi_an.A[i_cell][0] =  2*cell_x * pow(cell_y, 2.)
    grad_phi_an.A[i_cell][1] =  2*cell_y * pow(cell_x, 2.)

# check linear interpolation
phi_surf = fields.SurfaceScalarField(myMesh, "phi_surf")
for face_i in range(len(phi_surf.A)):
    phi_surf.A[face_i] = phi.A[face_i]
phi_surf_interp = calculated_fields.LinearFaceValue(phi)

# create grad_phi_num field with 2 non conjunctional iterations
grad_phi_num = calculated_fields.GaussCellGradient(phi, 2)

# write volume fields
volume_fields = numpy.array([phi, grad_phi_an, grad_phi_num])
output.write_mesh_file_with_fields(volume_fields)

# write surface fields
surface_fields = numpy.array([phi_surf, phi_surf_interp])
output.write_mesh_faces_with_field(surface_fields)

print ""
print "testGradient FINSHED"