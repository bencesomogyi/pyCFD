"""
2D test case for the Smith-Hutton problem
"""
__author__ = "Bence Somogyi"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Bence Somogyi"
__email__ = "bencesomogyi@ivt.tugraz.at"
__status__ = "Prototype"

# "global variables"
import pyCFD_config.config as Config_
Config_.__FIELDDIR__        = '_FIELD_FILES/'
Config_.__OUTDIR__          = '_OUTPUT/'
Config_.__OUTITERDIR__      = '_OUTPUT/ITERATIONS/'
Config_.__OUTITERDIRREL__   = 'ITERATIONS/'

import numpy
import pyCFD_mesh.readers as readers
import pyCFD_output.output as output
import pyCFD_fields.fields as fields
import pyCFD_fields.calculated_fields as calcfield
import pyCFD_calculation.time_loop as time_loop
import pyCFD_operators.generic_equation as generic_equation
#import pyCFD_operators.explicit_operators as explicit_operators
import pyCFD_operators.implicit_operators as implicit_operators

# create mesh object and read mesh data from MESHES directory
myMesh = readers.FoamMesh("smithhutton")

# create and initialize field U
U = fields.VectorField(myMesh, "U")

for i_cell in range(len(U.V)):
    cell_x = myMesh.cells[i_cell].C[0]
    cell_y = myMesh.cells[i_cell].C[1]
    U.V[i_cell][0] =  2. * cell_y * (1. - pow(cell_x,2))
    U.V[i_cell][1] = -2. * cell_x * (1. - pow(cell_y,2))

for patch_ in myMesh.patches:
    patch_length = len(patch_.faces)
    temp_vector = numpy.zeros((patch_length,3))
    for face_i,face_ in enumerate(patch_.faces):
        face_x = face_.C[0]
        face_y = face_.C[1]
        temp_vector[face_i][0] =  2. * face_y * (1. - pow(face_x,2))
        temp_vector[face_i][1] = -2. * face_x * (1. - pow(face_y,2))
    U.get_patch(patch_.name).set_patch_distributed(temp_vector, "fixedValue")
    
U.update_boundary_values()

# create rho field
rho = fields.ScalarField(myMesh, "rho", 1.0)

# create massflux field
m_dot = calcfield.MassFlux(U, rho)

# create scalar field with boundary conditions
phi = fields.ScalarField(myMesh, "phi", 1.0)
phi.get_patch("inlet1"      ).set_patch_uniform(1., "fixedValue"   )
phi.get_patch("outlet"      ).set_patch_uniform(0., "fixedGradient")
phi.get_patch("wall"        ).set_patch_uniform(0., "fixedGradient")
phi.get_patch("frontAndBack").set_patch_uniform(0., "fixedGradient")

# create time loop
save_fields = [U, phi]

start_time = 0.
stop_time  = 4.
time_step  = 0.05
save_step  = 1.

myLoop = time_loop.TimeLoop(save_fields, start_time, stop_time, time_step)
myLoop.uniform_save_times(save_step)

# set up equation for phi
phi_eqn = generic_equation.GenericScalarEquation(myMesh, phi, "bicg")

# save initial condition
if start_time == 0.:
    # clear output dirs
    output.clean_output_dirs()
    myLoop.save_current(0)

# iteration loop
for time_ in myLoop.times:
    myLoop.time = time_
    myLoop.print_time()
    
    phi_eqn.reset()
    
    phi_eqn += implicit_operators.DdtEuler(phi_eqn, myLoop.dt)    
#    phi_eqn += implicit_operators.Divergence(phi, m_dot, "UDS")
#    phi_eqn += implicit_operators.Divergence(phi, m_dot, "MINMOD")
    phi_eqn += implicit_operators.Divergence(phi, m_dot, "STOIC")
    
    phi_eqn.solve()
        
    myLoop.save_time()
    
phi_eqn.save_residual(myLoop.times,"phi")

print ""
print "testConvection_SmithHutton_UDS FINSHED"