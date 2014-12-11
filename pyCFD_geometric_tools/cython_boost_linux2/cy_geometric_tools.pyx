"""
cython module for calculations in the geomTools module
"""
__author__ = "Bence Somogyi"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Bence Somogyi"
__email__ = "bencesomogyi@ivt.tugraz.at"
__status__ = "Prototype"
# from __future__ import division
cimport numpy
import numpy
cimport cython

DTYPE = numpy.float
ctypedef numpy.float_t DTYPE_t

cdef extern from "math.h":
    float sqrt (float x)
    float abs (float x)

@cython.boundscheck(False)
# @cython.cdivision(True)

def cy_triangle_areas(numpy.ndarray[DTYPE_t, ndim=2] area_vector):
    """
    calculate triangle area from area vector
    
    :param area_vector: array of area vector
    :type area_vector:  float
    :return:            area
    :rtype:             float
    
    .. note:
        cython code to compile c library
    """
    cdef unsigned int n = len(area_vector)
    cdef numpy.ndarray[DTYPE_t, ndim=1] triangle_areas = numpy.zeros(n)
    for i in xrange(n):
        triangle_areas[i] = sqrt(area_vector[i,0]*area_vector[i,0] + area_vector[i,1]*area_vector[i,1] + area_vector[i,2]*area_vector[i,2])
    return triangle_areas

def cy_area_vector(numpy.ndarray[DTYPE_t, ndim=2] centroid, numpy.ndarray[DTYPE_t, ndim=2] coords_a, numpy.ndarray[DTYPE_t, ndim=2] coords_b):
    """
    calculate triangle area vector from coordinates of 3 vertices
    
    :param centroid: coordinate array of node 1
    :type centroid:  float
    :param coords_a: coordinate array of node 2
    :type coords_a:  float
    :param coords_b: coordinate array of node 3
    :type coords_b:  float
    :return:         array of area vector
    :rtype:          float
    
    .. note:
        cython code to compile c library
    """
    cdef unsigned int n = len(centroid)
    cdef numpy.ndarray[DTYPE_t, ndim=2] area_vector = numpy.zeros((n,3))
    cdef numpy.ndarray[DTYPE_t, ndim=1] vec_a = numpy.zeros(3)
    cdef numpy.ndarray[DTYPE_t, ndim=1] vec_b = numpy.zeros(3)
    for i in xrange(n):
        vec_a[0] = centroid[i,0] - coords_a[i,0]
        vec_a[1] = centroid[i,1] - coords_a[i,1]
        vec_a[2] = centroid[i,2] - coords_a[i,2]
        vec_b[0] = centroid[i,0] - coords_b[i,0]
        vec_b[1] = centroid[i,1] - coords_b[i,1]
        vec_b[2] = centroid[i,2] - coords_b[i,2]
        area_vector[i,0] = ((vec_a[1] * vec_b[2]) - (vec_a[2] * vec_b[1])) * 0.5
        area_vector[i,1] = ((vec_a[2] * vec_b[0]) - (vec_a[0] * vec_b[2])) * 0.5
        area_vector[i,2] = ((vec_a[0] * vec_b[1]) - (vec_a[1] * vec_b[0])) * 0.5
    return area_vector

def cy_value_close(DTYPE_t value, DTYPE_t compare):
    """
    cython implementation of numpy.allclose()
    
    Returns True if two values are equal within a tolerance.

    The tolerance value is positive, typically a very small number. The
    relative difference (rtol * abs(compare)) and the absolute difference atol
    are added together to compare against the absolute difference between value
    and compare.
    
    **Tolerances**
    
    * rtol = 1e-05
    
    * atol = 1e-08
    
    :param value:   value to be compared
    :type value:    float
    :param compare: value to compare with
    :type compare:  float
    
    .. note:
        cython code to compile c library
    """
    cdef DTYPE_t rtol = 1e-05
    cdef DTYPE_t atol = 1e-08
    if abs(value - compare) <= (atol + rtol * abs(compare)):
        return True
    return False

def cy_calc_cos(numpy.ndarray[DTYPE_t, ndim=1] areaVect, numpy.ndarray[DTYPE_t, ndim=1] vectElmFace):
    """
    calculate cosine between face are vector and the vector between the cell
    centroid and the face centroid to decide if area vector points inside or
    outside the cell
    
    :param areaVect:    coordinate array of the area vector
    :type areaVect:     float
    :param vectElmFace: coordinate array of the vector between cell and face centroids
    :type vectElmFace:  float
    :return:            cos value between the vectors
    :rtype:             float
    
    .. note:
        cython code to compile c library
    """
    cdef DTYPE_t areaVect_norm = sqrt(areaVect[0]*areaVect[0] + areaVect[1]*areaVect[1] + areaVect[2]*areaVect[2])
    cdef DTYPE_t vectElmFace_norm = sqrt(vectElmFace[0]*vectElmFace[0] + vectElmFace[1]*vectElmFace[1] + vectElmFace[2]*vectElmFace[2])
    cdef DTYPE_t cos_value = ( areaVect[0]*vectElmFace[0] + areaVect[1]*vectElmFace[1] + areaVect[2]*vectElmFace[2] ) / areaVect_norm / vectElmFace_norm
    return cos_value
