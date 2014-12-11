"""
module for calculating monitoring quantities
"""
__author__ = "Bence Somogyi"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Bence Somogyi"
__email__ = "bencesomogyi@ivt.tugraz.at"
__status__ = "Prototype"

import pyCFD_config.config as Config_
import numpy
import math

# "global variables"
__OUTITERDIR__    = Config_.__OUTITERDIR__
__OUTITERDIRREL__ = Config_.__OUTITERDIRREL__
__OUTDIR__        = Config_.__OUTDIR__

def cd(scalar_field, referece_velocity, rho, patch_name, flow_direction, mesh_width=1.):
    r"""
    calculate drag coefficient on a 2D mesh:
        
    :math:`c_d = \frac{\sum_f{\phi_f \left( \vec{A} \cdot \vec{e_{flow}} \right)}}{0.5 \rho v_{ref}^2}`    
    
    :param scalar_field:       field to monitor: :math:`\phi`
    :type scalar_field:        :class:`pyCFD_fields.fields.ScalarField`
    :param reference_velocity: reference velocity: :math:`v_{ref}`
    :type reference_velocity:  float
    :param rho:                reference density: :math:`\rho`
    :type rho:                 float
    :param patch_name:         patch on which integration should be performed
    :type patch_name:          string
    :param flow_direction:     array of reference flow direction: :math:`\vec{e_{flow}}`
    :type flow_direction:      numpy.array
    :param mesh_width:         default: 1.0, mesh width normal to the flow
    :type mesh_width:          float
    """
    if (scalar_field.type != "scalar"):
        raise TypeError("field should be a scalar field for calculating the cd")

    distributed_rho = True
    if isinstance(rho, float):
        distributed_rho = False
        
    field_patch = scalar_field.get_patch(patch_name)
    
    cd = 0.
    for face_ in field_patch.father[0].faces:
        cd_local = scalar_field.A[face_.id] * numpy.dot(face_.Sf, flow_direction)
        if distributed_rho:
            cd_local /= ( 0.5 * referece_velocity * referece_velocity * rho.A[face_.id] * mesh_width )
        cd += cd_local
    if distributed_rho == False:
        cd /= ( 0.5 * referece_velocity * referece_velocity * rho * mesh_width)
        
    return cd
    
def globalMass(massflux_field):
    mass_sum = 0.
    for face_ in massflux_field.father[0].faces:
        if face_.isBnd == False:
            continue
        mass_sum += massflux_field.A[face_.id]
        
    return mass_sum
    
def CFLMag(velocity_field, length_scale, dt):
    CFL_ = 0.
    for cell_ in velocity_field.father[0].cells:
        CFL_temp = math.sqrt(velocity_field.V[cell_.id][0]*velocity_field.V[cell_.id][0] \
                   + velocity_field.V[cell_.id][1]*velocity_field.V[cell_.id][1] \
                   + velocity_field.V[cell_.id][2]*velocity_field.V[cell_.id][2]) \
                   * dt / length_scale
        if CFL_temp > CFL_:
            CFL_ = CFL_temp
            
    return CFL_
    
def CFL(velocity_field, length_scale, dt):
    CFL_ = 0.
    for cell_ in velocity_field.father[0].cells:
        for comp_ in range(3):
            CFL_temp = abs(velocity_field.V[cell_.id][comp_]) * dt / length_scale
            if CFL_temp > CFL_:
                CFL_ = CFL_temp
            
    return CFL_
    