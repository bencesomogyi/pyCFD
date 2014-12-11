"""
2D flow around square cylinder
"""
__author__ = "Bence Somogyi"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Bence Somogyi"
__email__ = "bencesomogyi@ivt.tugraz.at"
__status__ = "Prototype"

import time
t00 = time.time()

# "global variables"
import pyCFD_config.config as Config_
Config_.__sparse__          = True
Config_.__FIELDDIR__        = '_FIELD_FILES/'
Config_.__OUTDIR__          = '_OUTPUT/'
Config_.__OUTITERDIR__      = '_OUTPUT/ITERATIONS/'
Config_.__OUTITERDIRREL__   = 'ITERATIONS/'
Config_.__IMEX__            = True

import numpy
import pyCFD_mesh.readers as readers
import pyCFD_output.output as output
import pyCFD_fields.fields as fields
import pyCFD_fields.calculated_fields as calcfield
import pyCFD_calculation.time_loop as time_loop
import pyCFD_operators.generic_equation as generic_equation
import pyCFD_operators.explicit_operators as explicit_operators
import pyCFD_operators.implicit_operators as implicit_operators
import pyCFD_monitors.monitors as monitors
import pyCFD_linear_solvers.linear_solvers as mylin
#import pyCFD_fields.initialization as init

# create mesh object and read mesh data from MESHES directory
myMesh = readers.FoamMesh("2dChannel_wide")
cell_volumes = myMesh.get_volumes()[:,0]
face_areas = myMesh.get_areas()

# create and initialize field U, V and W
U = fields.VectorField(myMesh, "U", numpy.array([1.,0.,0.]))
U.get_patch("inlet"       ).set_patch_uniform([1.,0.,0.], "fixedValue"   )
U.get_patch("outlet"      ).set_patch_uniform(0., "fixedGradient")
U.get_patch("frontAndRear").set_patch_uniform(0., "fixedGradient")
#U.load_init_fields("square_cyl_82p8s")
#U.load_init_fields()
U_old = fields.VectorField(myMesh, "U_old")
U_old.V = U.V
U_old.A = U.A

# create and initialize field p
p = fields.ScalarField(myMesh, "p", 0.)
p.get_patch("inlet"       ).set_patch_uniform(0., "fixedGradient")
p.get_patch("walls"       ).set_patch_uniform(0., "fixedGradient")
p.get_patch("frontAndRear").set_patch_uniform(0., "fixedGradient")
#p.load_init_fields("square_cyl_82p8s")
#p.load_init_fields()
# p_old for IMEX scheme
if Config_.__IMEX__:
    p_old = fields.ScalarField(myMesh, "p_old")
    p_old.V = p.V
    p_old.A = p.A

# create and initialize field p_corr
p_corr = fields.ScalarField(myMesh, "p_corr", 0.)
p_corr.get_patch("inlet"       ).set_patch_uniform(0., "fixedGradient")
p_corr.get_patch("walls"       ).set_patch_uniform(0., "fixedGradient")
p_corr.get_patch("frontAndRear").set_patch_uniform(0., "fixedGradient")
#p_corr.load_init_fields("square_cyl_82p8s")
#p_corr.load_init_fields()

# create rho field
rho = fields.ScalarField(myMesh, "rho", 1.)
one_over_rho = fields.ScalarField(myMesh, "1_rho", 1.)

# create massflux field
m_dot = calcfield.MassFlux(U, rho)
m_dot_old = fields.SurfaceScalarField(myMesh, "massFlux_old", )
m_dot_old.copy_field(m_dot)
#m_dot_mon = fields.SurfaceScalarField(myMesh, "massFlux_mon", )

# set Reynolds number
Re = 10.

# create time loop
div_field = fields.ScalarField(myMesh, "div", 0.)
save_fields = [U, p, p_corr, div_field]

# 1st run
start_time =   0.
stop_time  =   0.13
time_step  =   0.0065
save_step  =   0.065
vtk_start  =   0
 
myLoop = time_loop.TimeLoop(save_fields, start_time, stop_time, time_step)
myLoop.uniform_save_times(save_step)
myLoop.vtk_start = vtk_start

# set up equation for U and p_corr
U_eqn      = generic_equation.GenericVectorEquation(myMesh, U,       "gs"       , 0.7 )
p_corr_eqn = generic_equation.GenericScalarEquation(myMesh, p_corr,  "lu_solver", 1.  )
p_under_relaxation      = 1.

# save initial condition
if start_time == 0.:
    # clear output dirs
    output.clean_output_dirs()
    myLoop.save_current(0)

# prepare monitors
# linear solver residuals are added by default
# divergence
p_corr_eqn.add_monitor("divU_pred")
p_corr_eqn.add_monitor("divU_corr")
# drag coefficient
p_corr_eqn.add_monitor("cd")
if start_time == 0.:
    U_eqn.save_residual_header("U")
    p_corr_eqn.save_residual_header("p")
# global mass balance
p_corr_eqn.add_monitor("global_mass_cons")

# build Laplace term coefficient matrix now as it will not change
U_lapl = implicit_operators.LaplaceVec(U, 1./Re)
p_lapl = implicit_operators.Laplace(p_corr)

if start_time == 0.:
    print "\ncalculating lu decomposition for pressure..."
    tlu = time.time()
    if Config_.__sparse__:
        p_corr_eqn.p, p_corr_eqn.l, p_corr_eqn.u = mylin.lu_decomp(p_lapl.A.todense())
    else:
        p_corr_eqn.p, p_corr_eqn.l, p_corr_eqn.u = mylin.lu_decomp(p_lapl.A)
    print "DONE in "+str(time.time()-tlu)+" s"
    print "saving plu matrices..."
    p_corr_eqn.save_plu()
    print "DONE"
else:
    print "\nloading lu decomposition for pressure..."
    p_corr_eqn.load_plu("square_cyl_82p8s")
    print "DONE"

# iteration loop
for time_ in myLoop.times:
    myLoop.time = time_
    myLoop.dt = myLoop.find_variable_dt()
    myLoop.print_time()
    
    # update boundary conditions, reset matrices and right hand sides, update m_dot
    U_eqn.reset()
    p_corr_eqn.reset()
    p.update_boundary_values()
    
    if Config_.__IMEX__ == False:
        print "\nno momentum equation assembled!"
        input()
    
    if Config_.__IMEX__ == True:
        ## PREDICTOR IMEX
        ## convection, gradient: Adams-Bashfort
        ## diffusion:            Crank-Nicolson
        # assemble momentum equation
        U_eqn += implicit_operators.DdtEulerVec(U_eqn, myLoop.dt)
        Dt = fields.ScalarField(myMesh, "Dt_time")
        Dt.initialize_cell_with_vector(U_eqn.diag())
        U_eqn += explicit_operators.DivergenceVec(U,     m_dot,     "MINMOD") * 1.5
        U_eqn -= explicit_operators.DivergenceVec(U_old, m_dot_old, "MINMOD") * 0.5
        # use pre-calculated coefficient matrix for the laplace operator
        U_eqn -= implicit_operators.LaplaceVec(U, 1./Re, "", False) * 0.5
        U_eqn.A = U_eqn.A - U_lapl.A * 0.5
        U_eqn -= explicit_operators.LaplaceVec(U, 1./Re) * 0.5
        U_eqn += explicit_operators.Gradient(p,     one_over_rho) * 1.5
        U_eqn -= explicit_operators.Gradient(p_old, one_over_rho) * 0.5
    
    # solve momentum equation
    U_eqn.solve()
    
    # calculate predicted massflux field
    m_dot = calcfield.MassFlux(U, rho)
    
    ## CORRECTOR
    # correct m_dot with Rhie-Chow interpolation
    # skip this only after a restart!
    if time_ != myLoop.times[0] and time_ != 0.:
        D = fields.ScalarField(myMesh, "D")
        D.initialize_cell_with_vector(cell_volumes/U_eqn.diag())
        D_f = calcfield.LinearFaceValue(D)
        grad_p = calcfield.GaussCellGradient(p)
        grad_p_lin_f = calcfield.LinearFaceValue(grad_p)
        grad_p_f = calcfield.GaussFaceGradient(p)
        m_dot.A -= D_f.A * (grad_p_f.dot_Sf()[:,0] - grad_p_lin_f.dot_Sf()[:,0])
        m_dot_old_lin = calcfield.MassFlux(U_old, rho)
        Dt.V /= U_eqn.diag()
        Dt_f = calcfield.LinearFaceValue(Dt)
        m_dot.A += Dt_f.A * (m_dot_old.A - m_dot_old_lin.A)
    
    # monitor predictor divergence
    div = calcfield.Divergence(None, m_dot, "")
    print "predictor divergence: "+str(max(abs(div.V)))
    p_corr_eqn.append_to_monitor(max(abs(div.V)), "divU_pred")
    
    # assemble pressure-correction equation
    # use pre-calculated coefficient matrix for the laplace operator
    p_corr_eqn += implicit_operators.Laplace(p_corr, 1., "", False)
    p_corr_eqn.A = p_corr_eqn.A + p_lapl.A
    p_corr_eqn -= explicit_operators.Divergence(None, m_dot, "")
    p_corr_eqn.b[:,0] *= U_eqn.diag() / cell_volumes
    div_field.V[:] = p_corr_eqn.b[:,0]
    # solve pressure-correction equation
    p_corr_eqn.solve()
    
    # correct pressure
    p.V += p_corr.V * p_under_relaxation
    # correct velocity components
    p_corr.update_boundary_values()
    grad_p_corr = calcfield.GaussCellGradient(p_corr)
    U.V[:,0] -= grad_p_corr.V[:,0] * cell_volumes/U_eqn.diag() * p_under_relaxation
    U.V[:,1] -= grad_p_corr.V[:,1] * cell_volumes/U_eqn.diag() * p_under_relaxation
    U.V[:,2] -= grad_p_corr.V[:,2] * cell_volumes/U_eqn.diag() * p_under_relaxation
    U.update_boundary_values()
    # correct mass flux with Rhie-Chow interpolation
    if time_ != myLoop.times[0] and time_ != 0.:
        grad_p_corr_f = calcfield.GaussFaceGradient(p_corr)
        m_dot.A -= D_f.A * (grad_p_corr_f.dot_Sf()[:,0])
        m_dot_old.copy_field(m_dot)
    
    # overwrite old fields
    U_old.copy_field(U)
    if Config_.__IMEX__ == True:
        p_old = p
    
    # monitor remaining divergence
    div = calcfield.Divergence(None, m_dot, "")
    print "corrector divergence: "+str(max(abs(div.V)))
    p_corr_eqn.append_to_monitor(max(abs(div.V)), "divU_corr")
    
    # monitor drag coefficient
    current_cd = monitors.cd(p, 1., 1., "walls", numpy.array([1., 0., 0.]), 0.1)
    print "cd of cylinder: "+str(current_cd)
    p_corr_eqn.append_to_monitor(current_cd, "cd")
    
    # monitor global mass ballance
    mass_ballance = monitors.globalMass(m_dot)
    print "global mass ballance: "+str(mass_ballance)
    p_corr_eqn.append_to_monitor(mass_ballance, "global_mass_cons")
    
    # monitor CFL number
    CFL_ = monitors.CFL(U, 0.02, myLoop.dt)
    print "CFL: "+str(CFL_)
    
    # save results and residual data of time step
    myLoop.save_time()
    p_corr_eqn.append_current_residual(myLoop.time, "p")
    U_eqn.append_current_residual(myLoop.time, "U")
    
#    index_ = (abs(myLoop.saveTimes+myLoop.startTime-myLoop.time) < 1e-8).tolist().index(True) + 1 + myLoop.vtk_start    
#    output.write_mesh_faces_with_field(m_dot, str(index_)+"S")
    
# end of iteration loop

print "\nfinished square cylinder example in "+str(time.time()-t00)+" s"
