"""
module for field initialization methods
"""
__author__ = "Bence Somogyi"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Bence Somogyi"
__email__ = "bencesomogyi@ivt.tugraz.at"
__status__ = "Prototype"

import numpy
import sys
if sys.platform == 'win32':
    import pyCFD_operators.cython_boost_win32.cy_operators as cy_operators
elif sys.platform == 'linux2':
    import pyCFD_operators.cython_boost_linux2.cy_operators as cy_operators
else:
    sys.exit("unknown platform " + sys.platform + " found in geomTools.py, stopping...")

def init_linear_scalar_sphere_distribution(field, center_value, side_value,
                                           center_coords, side_distance):
    """
    Initialize a field with a spherical linear distribution. Cell and face
    centers are calculated based on distance from distribution center
        
    :param field: a scalar field
    :type field:  :class:`pyCFD_fields.fields.ScalarField`
    :param center_value: value at the center of the distribution
    :type center_value: float
    :param side_value: value at the side of the distribution
    :type side_value: float
    :param center_coords: coordinates of the center of the distribution
    :type center_coords: numpy.array
    :param side_distance: distance fro center to side
    :type side_distance: float
    """
    for cells_and_faces in range(2):
        #  treat cells first
        if cells_and_faces == 0:
            obj_list = field.father[0].cells
            do_cells = True
            do_faces = False
        #  than faces
        elif cells_and_faces == 1:
            obj_list = field.father[0].faces
            do_cells = False
            do_faces = True
        for obj_index, obj in enumerate(obj_list):
            obj_center = obj.C
#            distance = numpy.linalg.norm(numpy.add(center_coords, -obj_center))
            distance = cy_operators.cy_linalg_norm(numpy.add(center_coords, -obj_center))
            if distance < side_distance:
                current_value = center_value + (side_value - center_value) * distance / side_distance
            else:
                current_value = side_value
            if do_cells:
                field.V[obj_index] = current_value
            elif do_faces:
                field.A[obj_index] = current_value
                
def init_cell_values_in_box(field, box_min, box_max, init_value):
    """
    Re-initialize field values within a box selection
        
    :param field: a scalar field
    :type field:  :class:`pyCFD_fields.fields.ScalarField`
    :param box_min: lower bounding coordinates of the box
    :type box_min: numpy.array
    :param box_max: higher bounding coordinates of the box
    :type box_max: numpy.array
    :param init_value: desired initial value
    :type init_value: float
    """
    mesh_ = field.father[0]
    # select cells first
    selected_cells = []
    for cell_ in mesh_.cells:
        if cell_.C[0] > box_min[0] and cell_.C[0] < box_max[0] and cell_.C[1] > box_min[1] and cell_.C[1] < box_max[1] and cell_.C[2] > box_min[2] and cell_.C[2] < box_max[2]:
            selected_cells.append(cell_)
    # apply init value            
    for cell_ in selected_cells:
        field.V[cell_.id] = init_value