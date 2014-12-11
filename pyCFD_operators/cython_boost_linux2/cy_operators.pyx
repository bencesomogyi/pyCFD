# from __future__ import division
cimport numpy
import numpy
cimport cython

DTYPE = numpy.float
ctypedef numpy.float_t DTYPE_t

cdef extern from "math.h":
    float sqrt (float x)

@cython.boundscheck(False)
# @cython.cdivision(True)

def cy_linalg_norm(numpy.ndarray[DTYPE_t, ndim=1] vector):
    cdef unsigned int n = len(vector)
    cdef DTYPE_t norm = 0.0
    for i in xrange(n):
        norm = norm + vector[i]*vector[i]
    return sqrt(norm)
