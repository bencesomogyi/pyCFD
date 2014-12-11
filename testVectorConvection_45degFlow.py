"""
Convection test case, vector field transported with 45 deg flow
"""
__author__ = "Bence Somogyi"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Bence Somogyi"
__email__ = "bencesomogyi@ivt.tugraz.at"
__status__ = "Prototype"

import numpy
import pyCFD_mesh.readers as readers
import pyCFD_output.output as output
import pyCFD_fields.fields as fields
import pyCFD_fields.calculated_fields as calcfield
import pyCFD_calculation.time_loop as time_loop
import pyCFD_operators.generic_equation as generic_equation
import pyCFD_operators.explicit_operators as explicit_operators
import pyCFD_operators.implicit_operators as implicit_operators

# clear output dirs
output.clean_output_dirs()

# create mesh object and read mesh data from MESHES directory
myMesh = readers.FoamMesh()

# create and initialize field U
U = fields.VectorField(myMesh, "U")

for i_cell in range(len(U.V)):
    cell_x = myMesh.cells[i_cell].C[0]
    cell_y = myMesh.cells[i_cell].C[1]
    U.V[i_cell][0] = -1.
    U.V[i_cell][1] =  1.

for patch_ in myMesh.patches:
    patch_length = len(patch_.faces)
    temp_vector = numpy.zeros((patch_length,3))
    for face_i,face_ in enumerate(patch_.faces):
        face_x = face_.C[0]
        face_y = face_.C[1]
        temp_vector[face_i][0] = -1.
        temp_vector[face_i][1] =  1.
    U.get_patch(patch_.name).set_patch_distributed(temp_vector, "fixedValue")
U.update_boundary_values()

# create rho field
rho = fields.ScalarField(myMesh, "rho", 1.0)

# create massflux field
m_dot = calcfield.MassFlux(U, rho)

# create vector field with boundary conditions
phi = fields.VectorField(myMesh, "phi_vec", [1., 10., 100.])
phi.get_patch("phi0"        ).set_patch_uniform([0., 9., 90.], "fixedValue"   )
phi.get_patch("phi1"        ).set_patch_uniform([1., 10., 100.], "fixedValue"   )
phi.get_patch("outlet"      ).set_patch_uniform(0., "fixedGradient")
phi.get_patch("frontAndBack").set_patch_uniform(0., "fixedGradient")

# create time loop
save_fields = [U,phi]

start_time = 0.
stop_time  = 0.2
time_step  = 0.001
save_step  = 0.01

myLoop = time_loop.TimeLoop(save_fields, start_time, stop_time, time_step)
myLoop.uniform_save_times(save_step)

# set up equation for phi
phi_eqn = generic_equation.GenericVectorEquation(myMesh, phi, "bicg")

# save initial condition
myLoop.save_current(0)

# iteration loop
for time_ in myLoop.times:
    myLoop.time = time_
    myLoop.print_time()
    
    phi_eqn.reset()
    
    phi_eqn += implicit_operators.DdtEulerVec(phi, myLoop.dt)    
    phi_eqn += implicit_operators.DivergenceVec(phi, m_dot, "STOIC")
    
    phi_eqn.solve()
        
    myLoop.save_time()
    
phi_eqn.save_residual(myLoop.times)

print ""
print "testConvection_45degFlow_UDS FINSHED"