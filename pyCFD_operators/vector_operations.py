"""
module for fast vector operations
"""
__author__ = "Bence Somogyi"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Bence Somogyi"
__email__ = "bencesomogyi@ivt.tugraz.at"
__status__ = "Prototype"

import numpy

def VecCross(left, right):
    """
    array cross product re-implementation for vectors    
    """
    x = ((left[1] * right[2]) - (left[2] * right[1]))
    y = ((left[2] * right[0]) - (left[0] * right[2]))
    z = ((left[0] * right[1]) - (left[1] * right[0]))
    return numpy.array([x, y, z])
    
def VecDot(left, right):
    """
    array dot product re-implementation for vectors    
    """
    return left[0]*right[0] + left[1]*right[1] + left[2]*right[2]