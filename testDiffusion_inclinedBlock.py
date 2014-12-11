"""
2D test case diffusion in an inclined block
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

import pyCFD_mesh.readers as readers
import pyCFD_output.output as output
import pyCFD_fields.fields as fields
#import pyCFD_fields.calculated_fields as calcfield
import pyCFD_calculation.time_loop as time_loop
import pyCFD_operators.generic_equation as generic_equation
#import pyCFD_operators.explicit_operators as explicit_operators
import pyCFD_operators.implicit_operators as implicit_operators

# clear output dirs
output.clean_output_dirs()

# create mesh object and read mesh data from MESHES directory
myMesh = readers.FoamMesh("inclinedBlock")

# create scalar field with boundary conditions
phi = fields.ScalarField(myMesh, "phi", 273.)
phi.get_patch("phi0"        ).set_patch_uniform(273., "fixedValue"   )
phi.get_patch("phi1"        ).set_patch_uniform(300., "fixedValue"   )
phi.get_patch("frontAndRear").set_patch_uniform(0., "fixedGradient")

# create time loop
save_fields = [phi]

start_time = 0.
stop_time  = 3.
time_step  = 0.1
save_step  = 1.

myLoop = time_loop.TimeLoop(save_fields, start_time, stop_time, time_step)
myLoop.uniform_save_times(save_step)

# set up equation for phi
phi_eqn = generic_equation.GenericScalarEquation(myMesh, phi, "bicg")

# save initial condition
myLoop.save_current(0)

# iteration loop
for time_ in myLoop.times:
    myLoop.time = time_
    myLoop.print_time()
    
    phi_eqn.reset()
    
    phi_eqn += implicit_operators.DdtEuler(phi_eqn, myLoop.dt)    
    phi_eqn -= implicit_operators.Laplace(phi, 0.01)
    #phi_eqn -= implicit_operators.Laplace(phi, 0.01, "OVERRELAXED")
    
    phi_eqn.solve()
        
    myLoop.save_time()
    
phi_eqn.save_residual(myLoop.times)

print ""
print "testDiffusion_SmithHutton_explicit FINSHED"