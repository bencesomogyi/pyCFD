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

import pyCFD_mesh.readers as readers
import pyCFD_output.output as output
import pyCFD_fields.fields as fields
#import pyCFD_fields.calculated_fields as calcfield
import pyCFD_calculation.time_loop as time_loop
import pyCFD_operators.generic_equation as generic_equation
import pyCFD_operators.explicit_operators as explicit_operators
import pyCFD_operators.implicit_operators as implicit_operators


# clear output dirs
output.clean_output_dirs()

# create mesh object and read mesh data from MESHES directory
myMesh = readers.FoamMesh()

# create scalar field with boundary conditions
phi = fields.ScalarField(myMesh, "phi", 0.)
phi.get_patch("left"      ).set_patch_uniform(1., "fixedValue"      )
phi.get_patch("right"     ).set_patch_uniform(2., "fixedValue"      )
phi.get_patch("bottom"    ).set_patch_uniform(2., "fixedGradient"   )

# create time loop
save_fields = []
save_fields.append(phi)

start_time = 0.
stop_time  = 2.
time_step  = 1.
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
    
#    phi_eqn += implicit_operators.DdtEuler(phi, myLoop.dt)    
    phi_eqn -= implicit_operators.Laplace(phi)
    
    phi_eqn.solve()
        
    myLoop.save_time()
    
phi_eqn.save_residual(myLoop.times)

print ""
print "testDiffusion_SmithHutton_explicit FINSHED"