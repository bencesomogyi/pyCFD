"""
module for linear equation solvers
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
import scipy.sparse.linalg
import scipy.linalg
import scipy.io

## "global variables"
#import __builtin__
__sparse__ = Config_.__sparse__
#__sparse__ = True
if __sparse__:
    import scipy
    import scipy.sparse
    
if sys.platform == 'win32':
#    import pyCFD_linear_solvers.cython_boost_win32.cy_linear_solvers as cy_linear_solvers
    import cython_boost_win32.cy_linear_solvers as cy_linear_solvers
elif sys.platform == 'linux2':
    import pyCFD_linear_solvers.cython_boost_linux2.cy_linear_solvers as cy_linear_solvers
else:
    sys.exit("unknown platform " + sys.platform + " found in geomTools.py, stopping...")
  
def lu_decomp(A):
    """
    calculate LU decomposition of dense matrix A with pivoting
    
    :param A: array to be decomposed
    :type A:  numpy.array
    :return:  permutation matrix, lower triangular matrix, upper triangular matrix
    :rtype:   numpy.array
    """
    u         = numpy.array(A)#, dtype=float)
    n         = u.shape[1]
    l         = numpy.eye(n)
    p         = numpy.zeros((n,n))
    orig_rows = numpy.array(range(n))
    
    for j in range(0, n-1):
        # pivot
        new_index = abs(u[j:n,j]).argmax() + j
        if new_index != j:
            # exchange rows
            u[j,:], u[new_index,:] = numpy.array(u[new_index, :]), numpy.array(u[j,:])
            orig_rows[j], orig_rows[new_index] = orig_rows[new_index], orig_rows[j]
            if j != 0:
                l[j,0:j], l[new_index,0:j] = numpy.array(l[new_index,0:j]), numpy.array(l[j,0:j])
        for i in range(j+1, n):
            if u[j,j] == 0.:
                continue
            l[i,j] = u[i,j]/u[j,j]
            u[i,:] -= u[i,j]/u[j,j]*u[j,:]
            
    for j in range(n):
        p[orig_rows[j],j] = 1.
        
    return p,l,u
    
def lu_decomp_sparse(A):
    """
    calculate LU decomposition of sparse matrix A with pivoting
    
    :param A: array to be decomposed
    :type A:  scipy.sparse matrix
    :return:  permutation matrix, lower triangular matrix, upper triangular matrix
    :rtype:   scipy.sparse matrix
    """
    u         = scipy.sparse.dok_matrix(A)
    n         = u._shape[0]
    l         = scipy.sparse.eye(n,n)
    l         = scipy.sparse.dok_matrix(l)
    p         = scipy.sparse.dok_matrix((n, n))
    orig_rows = numpy.array(range(n))
    
    for j in range(0, n-1):
        # pivot
        u_nonzero = u.nonzero()
        filled_elements_in_column_j = (numpy.where(u_nonzero[1]==j))[0]
        ref = abs(u[j,j])
        new_index = j
        for i in filled_elements_in_column_j:
            i_ = u_nonzero[0][i]
            j_ = u_nonzero[1][i]
            if abs(u[i_,j_])>ref and i_>j:
                new_index = i_
                ref = u[i_,j_]
        if new_index != j:
            # exchange rows
            filled_elements_in_exchange_rows = numpy.append((numpy.where(u_nonzero[0]==j))[0], (numpy.where(u_nonzero[0]==new_index))[0])
            already_done = []
            j_list = [u_nonzero[1][k] for k in filled_elements_in_exchange_rows]
            # append alredy done u columns. As they have 0. elements, they are
            # not added to j_list, l matrix changes will be missing!
            j_list.extend([j__ for j__ in range(j)])
            for j_ in j_list:
                if j_ not in already_done:
                    already_done.append(j_)
                else:
                    continue
                if j_>=j:
                    u[j, j_], u[new_index, j_] = u[new_index, j_], u[j, j_]
                if j!=0 and j_<j:
                    l[j,j_], l[new_index,j_] = l[new_index,j_], l[j,j_]
            orig_rows[j], orig_rows[new_index] = orig_rows[new_index], orig_rows[j]
        u_nonzero = u.nonzero()
        filled_elements_in_column_j = (numpy.where(u_nonzero[1]==j))[0]
        i_list = [u_nonzero[0][i] for i in filled_elements_in_column_j]
        for i_ in i_list:
            if i_ <= j:
                continue
            if u[j,j] == 0.:
                continue
            if u[i_,j] == 0.:
                continue
            l[i_,j] = u[i_,j]/u[j,j]
            for k_ in range(j,n):
                u[i_,k_] -= l[i_,j]*u[j,k_]
            
    for j in range(n):
        p[orig_rows[j],j] = 1.
        
    return p,l,u

def lu_solver(A,b):
    """
    direct solver for linear system of equations using LU decomposition and
    backward and forward substituion.
    
    Calculation steps for the equation system :math:`A\phi=b`:
        
    * the coefficient matrix is decomposed: :math:`A=PLU`
        
    * forward substitution to solve for y: :math:`Ly=Pb`
    
    * backward substitution to solve for x: :math:`Ux=y`
    
    :param A: coefficient matrix
    :type A:  numpy.array
    :param b: right hand side
    :type b:  numpy.array
    :return:  solution vector of the equation system
    :rtype:   numpy.array
    """
    if __sparse__:
        p_,l_,u_ = lu_decomp_sparse(A)
        b_ = p_.transpose()
        b_ = b_.dot(b)
    else:
        p_,l_,u_ = lu_decomp(A)
        b_ = numpy.dot(p_.transpose(),b)
    n  = len(b)
    x  = numpy.zeros(n)
    y  =  numpy.zeros(n)
    for i in range(0,n):
        sum_ = 0.
        if i != 0:
            for j in range(0,i):
                sum_ += l_[i,j] * y[j]
        y[i] = (b_[i] - sum_) / l_[i,i]
        
    for i in range(n-1,-1,-1):
        sum_ = 0.
        if i != n-1:
            for j in range(n-1,i-1,-1):
                sum_ += u_[i,j] * x[j]
        x[i] = (y[i] - sum_) / u_[i,i]
        
    result = [x]
        
    return result
    

def lu_solver_plu(p_,l_,u_,b):
    """
    direct solver for linear system of equations using an existing LU
    decomposition of the coefficient matrix A. Result is obtained via
    backward and forward substituion. The substitution loops are implemented as
    cython functions in the
    :mod:`pyCFD_linear_solvers.cy_boost_linux2.cy_linear_solvers`:
    lu_solver_backward_loop and lu_solver_forward_loop.
    
    :param p_: permutation matrix
    :type p_:  numpy.array
    :param l_: lower triangular matrix
    :type l_:  numpy.array
    :param u_: upper triangular matrix
    :type u_:  numpy.array
    :param b:  right hand side
    :type b:   numpy.array
    :return:   solution vector of the equation system
    :rtype:    numpy.array
    """
    if __sparse__:
        b_ = p_.transpose()
        b_ = b_.dot(b)
    else:
        b_ = numpy.dot(p_.transpose(),b)
    b__ = b_.reshape(len(b_))
    n  = len(b)
    # x  = numpy.zeros(n)
    # y  =  numpy.zeros(n)
    # for i in range(0,n):
    #     sum = 0.
    #     if i != 0:
    #         for j in range(0,i):
    #             sum += l_[i,j] * y[j]
    #     y[i] = (b_[i] - sum) / l_[i,i]
    y = cy_linear_solvers.lu_solver_backward_loop(b__, l_, n)
        
    # for i in range(n-1,-1,-1):
    #     sum = 0.
    #     if i != n-1:
    #         for j in range(n-1,i-1,-1):
    #             sum += u_[i,j] * x[j]
    #     x[i] = (y[i] - sum) / u_[i,i]
    x = cy_linear_solvers.lu_solver_forward_loop(b__, y, u_, n)
        
    result = [x]
        
    return result

def gs(A,b,x0,tol,max_iter):
    """
    iterative solver for linear system of equations using Gauss-Seidel
    iterations. Two sub-iterations are performed in each iteration steps: once
    starting from front, once starting from rear.
    
    The sparse solution is implemented in cython.
    
    :param A:        coefficient matrix
    :type A:         numpy.array
    :param b:        right hand side
    :type b:         numpy.array
    :param x0:       initial condition vector
    :type x0:        numpy.array
    :param tol:      absolute tolerance between two iterations
    :type tol:       float
    :param max_iter: maximum number of iterations
    :type max_iter:  int
    :return:         solution vector of the equation system
    :rtype:          numpy.array
    """
    if __sparse__:
        A_nonzero = A.nonzero()
        row_i = A_nonzero[0]
        col_i = A_nonzero[1]
        mat_v = A.data
        dia_v = A.diagonal()
        b_ = b.reshape(len(b))
    n = A.shape[1]
    if x0 == None:
        x0 = numpy.zeros(n)
    new_x = 0.
    delta_x = 1.
    
    if __sparse__:
        x, iter_, delta_ = cy_linear_solvers.gs_sparse_loop(row_i,col_i,mat_v,dia_v,b_,x0,max_iter,tol,delta_x)
    else:
        iter_ = 0
        x = numpy.zeros(n)
        while iter_ < max_iter and tol < delta_x:
            for i in range(n):
                new_x = 1./A[i,i] * (b[i] - numpy.dot(A[i,:i], x[:i]) - numpy.dot(A[i,i+1:], x[i+1:]))
                if abs(new_x - x[i]) > delta_x:
                    delta_x = abs(new_x - x[i])
                x[i] = new_x
            for i in range(n-1,-1,-1):
                new_x = 1./A[i,i] * (b[i] - numpy.dot(A[i,:i], x[:i]) - numpy.dot(A[i,i+1:], x[i+1:]))
                if abs(new_x - x[i]) > delta_x:
                    delta_x = abs(new_x - x[i])
                x[i] = new_x
            iter_ += 1
    print "- gs solver reached "+str(delta_)+" after "+str(iter_)+" iterations"
    return [x]             
#    
#    
#    x = x0
#    iter = 0
#    while iter < max_iter and tol < delta_x:
#        if __sparse__:
#            for i in range(n):
#                sum_ = 0.
#                j_for_i = (numpy.where(A_nonzero[0]==i))[0]
#                for j_ in j_for_i:
#                    j = A_nonzero[1][j_]
#                    if j != i:
#                        sum_ += A[i,j]*x[j]
#                new_x = 1./A[i,i] * (b[i] - sum_)
#                if abs(new_x - x[i]) > delta_x:
#                    delta_x = abs(new_x - x[i])
#                x[i] = new_x                
#            for i in range(n-1,-1,-1):
#                sum_ = 0.
#                j_for_i = (numpy.where(A_nonzero[0]==i))[0]
#                for j_ in j_for_i:
#                    j = A_nonzero[1][j_]
#                    if j != i:
#                        sum_ += A[i,j]*x[j]
#                new_x = 1./A[i,i] * (b[i] - sum_)
#                if abs(new_x - x[i]) > delta_x:
#                    delta_x = abs(new_x - x[i])
#                x[i] = new_x
#        else:
#            for i in range(n):
#                new_x = 1./A[i,i] * (b[i] - numpy.dot(A[i,:i], x[:i]) - numpy.dot(A[i,i+1:], x[i+1:]))
#                if abs(new_x - x[i]) > delta_x:
#                    delta_x = abs(new_x - x[i])
#                x[i] = new_x
#            for i in range(n-1,-1,-1):
#                new_x = 1./A[i,i] * (b[i] - numpy.dot(A[i,:i], x[:i]) - numpy.dot(A[i,i+1:], x[i+1:]))
#                if abs(new_x - x[i]) > delta_x:
#                    delta_x = abs(new_x - x[i])
#                x[i] = new_x
#        iter += 1
#    return [x]
    

# lu test matrices
#A = scipy.io.mmread("../p_lapl_a.mtx")
#A = scipy.io.mmread("p_lapl_a.mtx")
#A = numpy.array([[2.,2.,0.,1.],[2.,0.,3.,2.],[4.,-3.,1.,1.],[6.,1.,-6.,0.]])
#A = numpy.array([[ 0.41667454,  0.54947262,  0.90061293,  0.77658974],
#                 [ 0.4386116,   0.5074603,   0.53885945,  0.45073554],
#                 [ 0.92089123,  0.57822875,  0.01764067,  0.54463878],
#                 [ 0.56181101,  0.24405833,  0.01401994,  0.25714378]])    
#A = numpy.random.rand(50,50)    
#A  = numpy.array([[3.,1.], [1.,2.]])
#A = numpy.load("../p_corr_eqn_A.npy")
#p_,l_,u_ = scipy.linalg.lu(A.todense())
#p_,l_,u_ = scipy.linalg.lu(A)
#A_ = A.tocsr()
#print A.todense()
#p_,l_,u_ = scipy.sparse.linalg.splu(A.todense())
#p,l,u = lu_decomp_sparse(A)
#print ""
#print "u:"
#print u.todense()
#print u_
#print (u_.todense()-u.todense()).max()
#print (u_.todense()-u.todense()).min()
#print ""
#print "diff l:"
#print l.todense()
#print l_
#print (l_-l.todense()).max()
#print (l_-l.todense()).min()
#print ""
#print "diff p:"
#print (p_-p.todense()).max()
#print (p_-p.todense()).min()
#print ""
#print numpy.dot(numpy.dot(p.todense(),l.todense()),u.todense())
#print ""
#print A
# test lu_solver
#b  = numpy.array([9.,8.])
#b = numpy.array([0.,-2.,-7.,6.])
#b = numpy.load("../p_corr_eqn_b.npy")
#x0 = numpy.linalg.solve(A, b)
#x1 = lu_solver(A,b)
#x2 = lu_solver_plu(p,l,u,b)
#for i_ in range(len(x0)):
#    print x0[i_], x2[i_]
#print x0
#print x1
#print x2
#print (x0-x1).max()
#print (x0-x2).max()

# test gs solver
#A = numpy.load("../U_eqn_A.npy")
#A_ = A
#A_ = scipy.sparse.dok_matrix(A)
#A_ = A_.tocsr()
#b = numpy.load("../U_eqn_bX.npy")
#A  = numpy.array([[3.,1.], [1.,2.]])
#A = numpy.array([[8.,2.,0.,1.],[2.,9.,3.,2.],[4.,-3.,5.,1.],[6.,1.,-6.,10.]])
#b = numpy.array([4.,-2.,-7.,6.])
#b  = numpy.array([9.,8.])
#x0 = numpy.linalg.solve(A,b)
#x1 = gs(A_,b,None,0.001,100)
#print x0
#print x1