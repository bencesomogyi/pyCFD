# from __future__ import division
import numpy
cimport numpy
cimport cython
from libc.math cimport abs

DTYPE = numpy.float
ctypedef numpy.float_t DTYPE_t

#cdef extern from "math.h":
#    float abs (float x)

@cython.boundscheck(False)
# @cython.cdivision(True)

def gs_sparse_loop(numpy.ndarray[int, ndim=1] row_indices,
                   numpy.ndarray[int, ndim=1] column_indices,
                   numpy.ndarray[DTYPE_t, ndim=1] matrix_values,
                   numpy.ndarray[DTYPE_t, ndim=1] diagonal_values,
                   numpy.ndarray[DTYPE_t, ndim=1] b,
                   numpy.ndarray[DTYPE_t, ndim=1] x0,
                   int max_iter,
                   DTYPE_t tol,
                   DTYPE_t delta_x                                ):
    """
    iterative solver for linear system of equations using Gauss-Seidel
    iterations. Two sub-iterations are performed in each iteration steps: once
    starting from front, once starting from rear.
    
    :param row_indices:     row indices of the non zero sparse matrix elements
    :type row_indices:      numpy.array
    :param column_indices:  column indices of the non zero sparse matrix elements
    :type column_indices:   numpy.array
    :param matrix_values:   non zero matrix coefficients
    :type matrix_values:    numpy.array
    :param diagonal_values: diagonal values of the sparse matrix
    :type diagonal_values:  numpy.array
    :param b:               right hand side of the equation system
    :type b:                numpy.array
    :param x0:              initial condition
    :type x0:               numpy.array
    :param max_iter:        maximum number of iterations
    :type max_iter:         int
    :return:                solution vector of the equation system
    :rtype:                 numpy.array
    """
    cdef unsigned int n = len(x0)
    cdef numpy.ndarray[DTYPE_t, ndim=1] x = x0
    cdef DTYPE_t sum_
    cdef DTYPE_t new_x = 0.
    cdef int iter_ = 0
    while iter_ < max_iter:
        delta_x = 0.
        for i in xrange(n):
            sum_ = 0.0
            j_for_i = (numpy.where(row_indices==i))[0]
            for j_ in j_for_i:
                j = column_indices[j_]
                if j != i:
                    sum_ = sum_ + matrix_values[j_]*x[j]
            new_x = 1./diagonal_values[i] * (b[i] - sum_)
            if abs(new_x - x[i]) > delta_x:
                delta_x = abs(new_x - x[i])
            x[i] = new_x
        if tol > delta_x:
            break
        delta_x = 0.
        for i in xrange(n-1,-1,-1):
            sum_ = 0.0
            j_for_i = (numpy.where(row_indices==i))[0]
            for j_ in j_for_i:
                j = column_indices[j_]
                if j != i:
                    sum_ = sum_ + matrix_values[j_]*x[j]
            new_x = 1./diagonal_values[i] * (b[i] - sum_)
            if abs(new_x - x[i]) > delta_x:
                delta_x = abs(new_x - x[i])
            x[i] = new_x
        if tol > delta_x:
            break
        iter_ += 1
    return x, iter_, delta_x

def lu_solver_backward_loop(numpy.ndarray[DTYPE_t, ndim=1] b_,
                            numpy.ndarray[DTYPE_t, ndim=2] l_,
                            int n                             ):  
    """
    backward substitution loop for the lu direct solver
    
    :param b_: right hand side
    :type b_:  numpy.array
    :param l_: lower triangular matrix
    :type l_:  numpy.array
    :param n:  length of the unknown vector
    :type n:   int
    :return:   intermediate solution
    :rtype:    numpy.array
    """  
    cdef numpy.ndarray[DTYPE_t, ndim=1] y = numpy.zeros(n, DTYPE)
    cdef DTYPE_t sum_
    for i in xrange(0,n):
        sum_ = 0.
        if i != 0:
            for j in range(0,i):
                sum_ += l_[i,j] * y[j]
        y[i] = (b_[i] - sum_) / l_[i,i]
    return y
    
def lu_solver_forward_loop(numpy.ndarray[DTYPE_t, ndim=1] b_,
                           numpy.ndarray[DTYPE_t, ndim=1] y,
                           numpy.ndarray[DTYPE_t, ndim=2] u_,
                           int n                             ):
    """
    forward substitution loop for the lu direct solver
    
    :param b_: right hand side
    :type b_:  numpy.array
    :param y:  intermediate solution
    :type y:   numpy.array
    :param n:  length of the unknown vector
    :type n:   int
    :return:   final solution
    :rtype:    numpy.array
    """  
    cdef numpy.ndarray[DTYPE_t, ndim=1] x = numpy.zeros(n, DTYPE)
    cdef DTYPE_t sum_
    for i in range(n-1,-1,-1):
        sum = 0.
        if i != n-1:
            for j in range(n-1,i-1,-1):
                sum += u_[i,j] * x[j]
        x[i] = (y[i] - sum) / u_[i,i]
    return x