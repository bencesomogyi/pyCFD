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
#import pyCFD_operators.explicit_operators as explicit_operators
import pyCFD_operators.implicit_operators as implicit_operators


# clear output dirs
output.clean_output_dirs()

# create mesh object and read mesh data from MESHES directory
myMesh = readers.FoamMesh()

# create vector field with boundary conditions
phi = fields.VectorField(myMesh, "phi_vec", [0., 9., 99.])
phi.get_patch("inlet2"      ).set_patch_uniform([1., 10., 100.], "fixedValue"   )
phi.get_patch("inlet1"      ).set_patch_uniform([0., 9., 99.], "fixedValue")
phi.get_patch("outlet"      ).set_patch_uniform([0., 9., 99.], "fixedValue")
phi.get_patch("wall"        ).set_patch_uniform([0., 9., 99.], "fixedValue")
phi.get_patch("frontAndBack").set_patch_uniform(0., "fixedGradient")

# create time loop
save_fields = []
save_fields.append(phi)

start_time = 0.
stop_time  = 1.
time_step  = 1.
save_step  = 1.

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
    
#    phi_eqn += implicit_operators.DdtEulerVec(phi, myLoop.dt)    
    phi_eqn -= implicit_operators.LaplaceVec(phi, 0.01)
    
    phi_eqn.solve()
        
    myLoop.save_time()
    
phi_eqn.save_residual(myLoop.times)

print ""
print "testVectorDiffusion_SmithHutton FINSHED"