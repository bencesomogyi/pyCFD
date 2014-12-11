"""
cython module for boosting general tasks
"""
__author__ = "Bence Somogyi"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Bence Somogyi"
__email__ = "bencesomogyi@ivt.tugraz.at"
__status__ = "Prototype"
# from __future__ import division
# cimport numpy
# import numpy
# cython: embedsignature=True
cimport cython

DTYPE = float
ctypedef float DTYPE_t

@cython.boundscheck(False)

def list_index(list myList, DTYPE_t value):
    """
    find index of defined value in a list
    
    :param myList: list of values
    :type myList:  float
    :param value:  value to search the index for
    :typevalue:    float
    
    .. note:
        cython code to compile c library
    """
    cdef unsigned int n = len(myList)
    for i in xrange(n):
        if myList[i] == value:
            return i
