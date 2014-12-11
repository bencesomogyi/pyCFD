"""
generic module for scalar and vector operators
"""
__author__ = "Bence Somogyi"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Bence Somogyi"
__email__ = "bencesomogyi@ivt.tugraz.at"
__status__ = "Prototype"

import pyCFD_config.config as Config_
import numpy

# "global variables"
#import __builtin__
__sparse__ = Config_.__sparse__
if __sparse__:
    import scipy

class GenericScalarOperator:
    """
    basic class for scalar operators
    """
    def __init__(self, mesh_):
        cell_number = len(mesh_.cells)
        if __sparse__:
            self.A = scipy.sparse.dok_matrix((cell_number, cell_number))
        else:
            self.A = numpy.zeros(cell_number)            
        """coefficient matrix of the linear equation system"""
        self.b = numpy.zeros((cell_number, 1))
        """vector of right hand side values"""
        self.father = mesh_
        """reference to mesh"""
        
    def __add__(self, other):
        """
        addition operator
        """
        if __sparse__:
            self.A = self.A + other.A
        else:
            self.A = numpy.add(self.A, other.A)
        self.b = numpy.add(self.b, other.b)
        return self
        
    def __sub__(self, other):
        """
        substraction operator
        """
        if __sparse__:
            self.A = self.A - other.A
        else:
            self.A = numpy.subtract(self.A, other.A)
        self.b = numpy.subtract(self.b, other.b)
        return self
        
    def __neg__(self):
        """
        return negated
        """
        if __sparse__:
            self.A = -self.A
        else:
            self.A = numpy.negative(self.A)
        self.b = numpy.negative(self.b)
        return self
        
    def __mul__(self, other):
        """
        multiplication operator
        """
        if isinstance(other, float):
            self.A = self.A * other
            self.b = self.b * other
            return self
        else:
            raise TypeError(
            "unsupported operand type(s) for *: '{}' and '{}'"
            ).format(self.__class__, type(other))
            
    def __div__(self, other):
        """
        division by float
        """
        if isinstance(other, float):
            self.A = self.A / other
            self.b = self.b / other
            return self
        else:
            raise TypeError(
            "unsupported operand type(s) for /: '{}' and '{}'"
            ).format(self.__class__, type(other))
         
    def __rmul__(self, other):
        """
        right multiplication for floats
        """
        if isinstance(other, float):
            result   = GenericScalarOperator(self.father)
            result.A = self.A * other
            result.b = self.b * other
            return result
        else:
            raise TypeError(
            "unsupported operand type(s) for left*: '{}' and '{}'"
            ).format(self.__class__, type(other))
            
    def __rdiv__(self, other):
        """
        right division for floats
        """
        if isinstance(other, float):
            self.A = numpy.divide(other*numpy.ones(self.A.shape), self.A)
            self.b = numpy.divide(other*numpy.ones(self.b.shape), self.b)
            return self
        else:
            raise TypeError(
            "unsupported operand type(s) for left/: '{}' and '{}'"
            ).format(self.__class__, type(other))
        
    def fix_cell_value(self, cell_index):
        """
        modify the equations coefficient matrix to fix the value in a cell.
        Coefficient in the main diagonal is set to 1.0 and other coeffiecient
        in the cells row are set to 0.0.
        
        :param cell_index: index of the cell where value should be fixed
        :type cell_index:  int
        """
        cell_number = len(self.father.cells)
        for i_ in xrange(cell_number):
            self.A[cell_index, i_] = 0.
        self.A[cell_index, cell_index] = 1.
            
class GenericVectorOperator:
    """
    basic class for vector operators
    """
    def __init__(self, mesh_):
        cell_number = len(mesh_.cells)
        if __sparse__:
            self.A  = scipy.sparse.dok_matrix((cell_number, cell_number))
        else:
            self.A = numpy.zeros(cell_number)
        """coefficient matrix of the linear equation system for the x component"""
        self.bX = numpy.zeros((cell_number, 1))
        """vector of right hand side values for the x component"""
        self.bY = numpy.zeros((cell_number, 1))
        """vector of right hand side values for the y component"""
        self.bZ = numpy.zeros((cell_number, 1))
        """vector of right hand side values for the z component"""
        self.father = mesh_
        """reference to mesh"""
        
    def __add__(self, other):
        """
        addition operator
        """
        if __sparse__:
            self.A  = self.A + other.A
        else:
            self.A  = numpy.add(self.A, other.A)
        self.bX = numpy.add(self.bX, other.bX)
        self.bY = numpy.add(self.bY, other.bY)
        self.bZ = numpy.add(self.bZ, other.bZ)
        return self
        
    def __sub__(self, other):
        """
        substraction operator
        """
        self.A  = numpy.subtract(self.A, other.A)
        self.bX = numpy.subtract(self.bX, other.bX)
        self.bY = numpy.subtract(self.bY, other.bY)
        self.bZ = numpy.subtract(self.bZ, other.bZ)
        return self
        
    def __neg__(self):
        """
        return negated
        """
        self.A  = numpy.negative(self.A)
        self.bX = numpy.negative(self.bX)
        self.bY = numpy.negative(self.bY)
        self.bZ = numpy.negative(self.bZ)
        return self
        
    def __mul__(self, other):
        """
        multiplication operator
        """
        if isinstance(other, float):
            self.A  *= other
            self.bX *= other
            self.bY *= other
            self.bZ *= other
            return self
        else:
            raise TypeError(
            "unsupported operand type(s) for *: '{}' and '{}'"
            ).format(self.__class__, type(other))
            
    def __div__(self, other):
        """
        division by float
        """
        if isinstance(other, float):
            self.A  = self.A / other
            self.bX = self.bX / other
            self.bY = self.bY / other
            self.bZ = self.bZ / other
            return self
        else:
            raise TypeError(
            "unsupported operand type(s) for /: '{}' and '{}'"
            ).format(self.__class__, type(other))
         
    def __rmul__(self, other):
        """
        right multiplication for floats
        """
        if isinstance(other, float):
            result    = GenericVectorOperator(self.father)
            result.A  = self.A * other
            result.bX = self.bX * other
            result.bY = self.bY * other
            result.bZ = self.bZ * other
            return result
        else:
            raise TypeError(
            "unsupported operand type(s) for left*: '{}' and '{}'"
            ).format(self.__class__, type(other))
            
    def __rdiv__(self, other):
        """
        right division for floats
        """
        if isinstance(other, float):
            self.A  = numpy.divide(other*numpy.ones(self.A.shape), self.A)
            self.bX = numpy.divide(other*numpy.ones(self.bX.shape), self.bX)
            self.bY = numpy.divide(other*numpy.ones(self.bY.shape), self.bY)
            self.bZ = numpy.divide(other*numpy.ones(self.bZ.shape), self.bZ)
            return self
        else:
            raise TypeError(
            "unsupported operand type(s) for left/: '{}' and '{}'"
            ).format(self.__class__, type(other))                       