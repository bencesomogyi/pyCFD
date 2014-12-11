"""
module for generic scalar and vector equations
"""
__author__ = "Bence Somogyi"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Bence Somogyi"
__email__ = "bencesomogyi@ivt.tugraz.at"
__status__ = "Prototype"

import pyCFD_config.config as Config_
import sys
import numpy
import scipy.sparse.linalg as lin
import pyCFD_operators.generic_operator
import pyCFD_output.output
import pyCFD_linear_solvers.linear_solvers as mylin

# "global variables"
#import __builtin__
__OUTDIR__ = Config_.__OUTDIR__
__sparse__ = Config_.__sparse__
if __sparse__:
    import scipy

class GenericScalarEquation(pyCFD_operators.generic_operator.GenericScalarOperator):
    """
    basic class for scalar equations
    
    Available solver types:
        
    * conjugate gradient: :class:`numpy.linalg.cg` - python built-in
    
    * biconjugate gradient: :class:`numpy.linalg.bicg` - python built-in
    
    * Gauss-Seidel: :class:`pyCFD_linear_solvers.linear_solvers.gs` - own implementation
    
    * LU solver: :class:`pyCFD_linear_solvers.linear_solvers.lu_solver_plu` - own implementation
    """
    def __init__(self, mesh_, volume_field, solver_type, under_relax = 1.):
        """
        **Constructor**
        
        :param mesh_:        mesh object
        :type mesh_:         :class:`pyCFD_mesh.generic_mesh.GenericMesh`
        :param volume_field: variable scalar field
        :type volume_field:  :class:`pyCFD_fields.fields.ScalarField`
        :param solver_type:  solver type to be used
        :type solver_type:   string
        :param under_relax:  default: 1.0, explicit underrelaxation factor
        :type under_relax:   float
        """
        pyCFD_operators.generic_operator.GenericScalarOperator.__init__(self, mesh_)
        self.x = []
        """vector of unknowns"""
        self.x_old = volume_field.V
        """field vector from previsus iteration/time step"""
        self.x_old_old = volume_field.V
        """field vector from previsus previous iteration/time step"""
        self.solver = solver_type
        """string with iterative solver type"""
        self.field = volume_field
        """reference to field"""
        self.residuals = []
        """list of residual values"""        
        self.monitors = []
        """list of monitor values"""        
        self.monitor_names = []
        """list of monitor names"""
        self.times = []
        """list of time steps"""
        self.under_relax = under_relax
        """under relaxation for the solution"""
        self.p = []
        """permutation matrix for lu solver"""
        self.l = []
        """lower triangluation matrix for lu solver"""
        self.u = []
        """upper triangluation matrix for lu solver"""
        
        # check solver type
        available_solvers = ["cg", "bicg", "lu_solver", "gs"]
        if solver_type not in available_solvers:
            print "not supported linear solver "+solver_type+" defined, availabe types are:"
            print available_solvers
            sys.exit()
            
    def solve(self, x_0=None, tol=1e-05, max_iter=100, exchange_zero=1e-16):
        """
        Solve the equation *A * x = b* using an linear equation solver. After
        solving old unknown vectors and volume field values are overwritten.
        """
        if self.solver   == "cg":
            solver_return = lin.cg(self.A, self.b, x_0, tol, max_iter)
        elif self.solver == "bicg":
            solver_return = lin.bicg(self.A, self.b, x_0, tol, max_iter)
        elif self.solver == "lu_solver":
             solver_return = mylin.lu_solver_plu(self.p, self.l, self.u, self.b)
        elif self.solver == "lu_solver":
             solver_return = mylin.gs(self.A, self.b, x_0, tol, max_iter)
        self.x         = solver_return[0]
        self.x_old_old = self.x_old
        
        # find max absolute residual
        diff_ = abs(self.x - self.x_old)
        max_diff_ = max(diff_)
        if max_diff_ == 0.:
            self.residuals.append(exchange_zero)
        else:
            self.residuals.append(max_diff_)
        print 'Residual for ' + self.field.name + ': ' + str(max(diff_))
        
        #update volume field
        self.field.V = self.x_old + self.under_relax*(self.x - self.x_old)
#        self.field.V = self.x
        
        return True
        
    def save_residual(self, times, name=""):
        """
        save all stored residuals and other monitoring quantities to
        :mod:`pyCFD_config.config`.__OUTDIR__ at the end of the calculation as CSV file
        
        saved file is residuals_<field_name>.csv
        """
        if self.monitors != []:
            headers_ = []
            headers_.append('time')
            headers_.append('res_'+self.field.name)
            vectors_ = []
            vectors_.append(times)
            vectors_.append(self.residuals)
            for i_,mon_ in enumerate(self.monitor_names):
                headers_.append(mon_)
                vectors_.append(self.monitors[i_])
            pyCFD_output.output.write_csv_file(headers_, vectors_, __OUTDIR__+'residuals_'+name+'.csv')
        else:
            pyCFD_output.output.write_csv_file(['time', 'res_'+self.field.name], [times, self.residuals], __OUTDIR__+'residuals_'+name+'.csv')
            
    def save_residual_header(self, name=""):
        """
        save only headers of residuals and monitoring values to
        :mod:`pyCFD_config.config`.__OUTDIR__ as CSV file
        
        saved file is residuals_<field_name>.csv
        """
        headers_ = ['time', 'res_'+self.field.name]
        if self.monitors != []:
            for mon_ in self.monitor_names:
                headers_.append(mon_)
        pyCFD_output.output.write_row_to_csv_file(headers_, __OUTDIR__+'residuals_'+name+'.csv')
                
    def append_current_residual(self, time, name=""):
        """
        append current residuals and monitoring values to the CSV file
        
        saved file is residuals_<field_name>.csv
        """
        vectors_ = [time, self.residuals[-1]]
        if self.monitors != []:
            for i_ in xrange(len(self.monitor_names)):
                vectors_.append(self.monitors[i_][-1])
        pyCFD_output.output.write_row_to_csv_file(vectors_, __OUTDIR__+'residuals_'+name+'.csv', 'a')
        
    def relax(self):
        """
        Relax the solution of solve().
        
        .. warning:
            This function is not implemented.
        """
        return True
        
    def reset(self):
        """
        Reset coefficient matrix and RHS, update boundary face values and
        update x_old with current field.
        """
        if __sparse__:
            self.A = scipy.sparse.dok_matrix(self.A.shape)
        else:
            self.A = numpy.zeros(self.A.shape)
        self.b = numpy.zeros(self.b.shape)
        self.field.update_boundary_values()
        self.x_old = self.field.V
            
    def write_A(self, name_=""):
        """
        save the equations coefficient matrix to :mod:`pyCFD_config.config`.__OUTDIR__ as DAT file.
        """
        if name_ == "":
            name_ = __OUTDIR__ + self.field.name + '_A.dat'
        numpy.savetxt(name_, self.A)
            
    def write_b(self, name_=""):
        """
        save the equations right hand side to :mod:`pyCFD_config.config`.__OUTDIR__ as DAT file.
        """
        if name_ == "":
            name_ = __OUTDIR__ + self.field.name + '_b.dat'
        numpy.savetxt(name_, self.b)
        
    def add_monitor(self, name_):
        """
        add a new monitoring quantity to the residuals
        
        :param name_: name of new monitoring value
        :type name_:  string
        """
        self.monitors.append([])
        self.monitor_names.append(name_)
        
    def append_to_monitor(self, value_, name_):
        """
        append new value to an existing monitoring quantity
        
        :param: value_: new value
        :type value_:   float
        :param name_:   name of existing monitoring quantity
        :type name_:    string
        """
        i_ = self.monitor_names.index(name_)
        self.monitors[i_].append(value_)
        
    def fix_cell_value(self, cell_index):
        """
        modify the equations coefficient matrix to fix the value in a cell.
        Coefficient in the main diagonal is set to 1.0 and other coeffiecient
        in the cells row are set to 0.0.
        
        :param cell_index: index of the cell where value should be fixed
        :type cell_index:  int
        """
        cell_number = len(self.field.mesh.cells)
        for i_ in xrange(cell_number):
            self.A[cell_index, i_] = 0.
        self.A[cell_index, cell_index] = 1.
        
    def save_plu(self):
        """
        save p, l and u matrices to :mod:`pyCFD_config.config`.__FIELDDIR__
        """
        numpy.save(Config_.__FIELDDIR__+self.field.name+"_p", self.p)
        numpy.save(Config_.__FIELDDIR__+self.field.name+"_l", self.l)
        numpy.save(Config_.__FIELDDIR__+self.field.name+"_u", self.u)
        
    def load_plu(self, dir_name=""):
        """
        save p, l and u matrices to :mod:`pyCFD_config.config`.__FIELDDIR__
        """
        if dir_name == "":
            self.p = numpy.load(Config_.__FIELDDIR__+self.field.name+"_p.npy")
            self.l = numpy.load(Config_.__FIELDDIR__+self.field.name+"_l.npy")
            self.u = numpy.load(Config_.__FIELDDIR__+self.field.name+"_u.npy")
        else:
            self.p = numpy.load(Config_.__FIELDDIR__+"/"+dir_name+"/"+self.field.name+"_p.npy")
            self.l = numpy.load(Config_.__FIELDDIR__+"/"+dir_name+"/"+self.field.name+"_l.npy")
            self.u = numpy.load(Config_.__FIELDDIR__+"/"+dir_name+"/"+self.field.name+"_u.npy")
        
    def diag(self):
        """
        return matrix main diagonal as vector
        """
        mesh_ = self.father
        cell_number = len(mesh_.cells)
        diagonal_ = numpy.array([self.A[i_, i_] for i_ in xrange(cell_number)])
        return diagonal_
        
class GenericVectorEquation(pyCFD_operators.generic_operator.GenericVectorOperator):
    """
    basic class for vector equations
    
    Available solver types:
        
    * conjugate gradient: :class:`numpy.linalg.cg` - python built-in
    
    * biconjugate gradient: :class:`numpy.linalg.bicg` - python built-in
    
    * Gauss-Seidel: :class:`pyCFD_linear_solvers.linear_solvers.gs` - own implementation
    
    * LU solver: :class:`pyCFD_linear_solvers.linear_solvers.lu_solver_plu` - own implementation
    """
    def __init__(self, mesh_, volume_field, solver_type, under_relax = 1.):
        """
        **Constructor**
        
        :param mesh_:        mesh object
        :type mesh_:         :class:`pyCFD_mesh.generic_mesh.GenericMesh`
        :param volume_field: variable scalar field
        :type volume_field:  :class:`pyCFD_fields.fields.ScalarField`
        :param solver_type:  solver type to be used
        :type solver_type:   string
        :param under_relax:  default: 1.0, explicit underrelaxation factor
        :type under_relax:   float
        """
        pyCFD_operators.generic_operator.GenericVectorOperator.__init__(self, mesh_)
        self.xX = []
        """vector of unknowns of x component"""
        self.xY = []
        """vector of unknowns of y component"""
        self.xZ = []
        """vector of unknowns of z component"""
        self.xX_old = volume_field.V[:,0]
        """field vector of x component from previsus iteration/time step"""
        self.xY_old = volume_field.V[:,1]
        """field vector of y component from previsus iteration/time step"""
        self.xZ_old = volume_field.V[:,2]
        """field vector of z component from previsus iteration/time step"""
        self.xX_old_old = volume_field.V[:,0]
        """field vector of x component from previsus previous iteration/time step"""
        self.xY_old_old = volume_field.V[:,1]
        """field vector of y component from previsus previous iteration/time step"""
        self.xZ_old_old = volume_field.V[:,2]
        """field vector of z component from previsus previous iteration/time step"""
        self.solver = solver_type
        """string with iterative solver type - "cg"/"bicg" """
        self.field = volume_field
        """reference to field"""
        self.residualsX = []
        """list of residual values for x component"""
        self.residualsY = []
        """list of residual values for x component"""
        self.residualsZ = []
        """list of residual values for x component"""
        self.times = []
        """list of time steps"""
        self.under_relax = under_relax
        """under relaxation for the solution"""
        
        # check solver type
        available_solvers = ["cg", "bicg", "gs"]
        if solver_type not in available_solvers:
            print "not supported linear solver "+solver_type+" defined, availabe types are:"
            print available_solvers
            sys.exit()
            
    def solve(self, x_0=None, tol=1e-5, max_iter=100, exchange_zero=1e-16):
        """
        Solve the equation *A * x = b* using an iterative solver. After solving
        old unknown vectors and volume field values are overwritten.
        """
        if x_0 == None:
            x_0X = x_0
            x_0Y = x_0
            x_0Z = x_0
        else:
            x_0X = x_0[:,0]
            x_0Y = x_0[:,1]
            x_0Z = x_0[:,2]
        if self.solver == "cg":
            solver_return_x = lin.cg(self.A, self.bX, x_0X, tol, max_iter)
            solver_return_y = lin.cg(self.A, self.bY, x_0Y, tol, max_iter)
            solver_return_z = lin.cg(self.A, self.bZ, x_0Z, tol, max_iter)
        elif self.solver == "bicg":
            solver_return_x = lin.bicg(self.A, self.bX, x_0X, tol, max_iter)
            solver_return_y = lin.bicg(self.A, self.bY, x_0Y, tol, max_iter)
            solver_return_z = lin.bicg(self.A, self.bZ, x_0Z, tol, max_iter)
        elif self.solver == "gs":
            solver_return_x = mylin.gs(self.A, self.bX, x_0X, tol, max_iter)
            solver_return_y = mylin.gs(self.A, self.bY, x_0Y, tol, max_iter)
            solver_return_z = mylin.gs(self.A, self.bZ, x_0Z, tol, max_iter)
        
        self.xX         = solver_return_x[0]
        self.xX_old_old = self.xX_old
        
        self.xY         = solver_return_y[0]
        self.xY_old_old = self.xY_old
        
        self.xZ         = solver_return_z[0]
        self.xZ_old_old = self.xZ_old
        
        # find max absolute residual for x
        diff_x = abs(self.xX - self.xX_old)
        max_diff_x = max(diff_x)
        if max_diff_x == 0.:
            self.residualsX.append(exchange_zero)
        else:
            self.residualsX.append(max_diff_x)
        print 'Residual for ' + self.field.name + ' x: ' + str(max(diff_x))
        
        # find max absolute residual for y
        diff_y = abs(self.xY - self.xY_old)
        max_diff_y = max(diff_y)
        if max_diff_y == 0.:
            self.residualsY.append(exchange_zero)
        else:
            self.residualsY.append(max_diff_y)
        print 'Residual for ' + self.field.name + ' y: ' + str(max(diff_y))
        
        # find max absolute residual for z
        diff_z = abs(self.xZ - self.xZ_old)
        max_diff_z = max(diff_z)
        if max_diff_z == 0.:
            self.residualsZ.append(exchange_zero)
        else:
            self.residualsZ.append(max_diff_z)
        print 'Residual for ' + self.field.name + ' z: ' + str(max(diff_z))
        
        #update volume field
        self.field.V[:,0] = self.xX_old + self.under_relax*(self.xX - self.xX_old)
        self.field.V[:,1] = self.xY_old + self.under_relax*(self.xY - self.xY_old)
        self.field.V[:,2] = self.xZ_old + self.under_relax*(self.xZ - self.xZ_old)
        
        return True
        
    def save_residual(self, times, name=""):
        """
        save all stored residuals and other monitoring quantities to
        :mod:`pyCFD_config.config`.__OUTDIR__ at the end of the calculation as CSV file
        
        saved file is residuals_<field_name>.csv
        """
        pyCFD_output.output.write_csv_file(['time', 'res_'+self.field.name+'_X', 'res_'+self.field.name+'_Y', 'res_'+self.field.name+'_Z'], [times, self.residualsX, self.residualsY, self.residualsZ], __OUTDIR__+'residuals_'+name+'.csv')
                    
    def save_residual_header(self, name=""):
        """
        save only headers of residuals and monitoring values to
        :mod:`pyCFD_config.config`.__OUTDIR__ as CSV file
        
        saved file is residuals_<field_name>.csv
        """
        headers_ = ['time', 'res_'+self.field.name+'_X', 'res_'+self.field.name+'_Y', 'res_'+self.field.name+'_Z']
        pyCFD_output.output.write_row_to_csv_file(headers_, __OUTDIR__+'residuals_'+name+'.csv')
                
    def append_current_residual(self, time, name=""):
        """
        append current residuals and monitoring values to the CSV file
        
        saved file is residuals_<field_name>.csv
        """
        vectors_ = [time, self.residualsX[-1], self.residualsY[-1], self.residualsZ[-1]]
        pyCFD_output.output.write_row_to_csv_file(vectors_, __OUTDIR__+'residuals_'+name+'.csv', 'a')
        
    def reset(self):
        """
        Reset coefficient matrix and RHS, update boundary face values and
        update x_old with current field.
        """
        if __sparse__:
            self.A = scipy.sparse.dok_matrix(self.A.shape)
        else:
            self.A = numpy.zeros(self.A.shape)
        self.bX = numpy.zeros(self.bX.shape)
        self.bY = numpy.zeros(self.bY.shape)
        self.bZ = numpy.zeros(self.bZ.shape)
        self.field.update_boundary_values()
        self.xX_old = self.field.V[:,0]
        self.xY_old = self.field.V[:,1]
        self.xZ_old = self.field.V[:,2]
        
    def write_A(self, name_=""):
        """
        save the equations coefficient matrix to :mod:`pyCFD_config.config`.__OUTDIR__ as DAT file.
        """
        if name_ == "":
            name_ = __OUTDIR__ + self.field.name + '_A.dat'
        numpy.savetxt(name_, self.A)
            
    def write_b(self, name_=""):
        """
        save the component equations right hand sides to :mod:`pyCFD_config.config`.__OUTDIR__ as
        DAT files.
        """
        if name_ == "":
            name_x = __OUTDIR__ + self.field.name + '_bX.dat'
            name_y = __OUTDIR__ + self.field.name + '_bY.dat'
            name_z = __OUTDIR__ + self.field.name + '_bZ.dat'
        numpy.savetxt(name_x, self.bX)
        numpy.savetxt(name_y, self.bY)
        numpy.savetxt(name_z, self.bZ)
        
    def diag(self):
        """
        return matrix main diagonal as vector
        """
        mesh_ = self.father
        cell_number = len(mesh_.cells)
        diagonal_ = numpy.array([self.A[i_, i_] for i_ in xrange(cell_number)])
        return diagonal_