"""
module for explicit operators
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
import pyCFD_fields.calculated_fields
import pyCFD_operators.generic_operator
if sys.platform == 'win32':
    import pyCFD_operators.cython_boost_win32.cy_operators as cy_operators
elif sys.platform == 'linux2':
    import pyCFD_operators.cython_boost_linux2.cy_operators as cy_operators
else:
    sys.exit("unknown platform " + sys.platform + " found in geomTools.py, stopping...")

# "global variables"
#import __builtin__
__sparse__ = Config_.__sparse__
if __sparse__:
    import scipy
          
class Divergence(pyCFD_operators.generic_operator.GenericScalarOperator):
    r"""
    Explicit divergence operator for :class:`pyCFD_fields.fields.ScalarField`
    using Picard iteration with known mass fluxes. Divergence is calculated
    using the Gauss theorem:
    
    .. math::
        
        \int_{V} \nabla \left( \rho \vec{v} \phi \right) dV &= \oint_{A} \left( \rho \vec{v} \phi \right) \cdot \vec{dA} \\
                                                               &= \sum_{f} \dot{m}_f * \phi_f
                                                               
    The explicit divergence operator returns the explicit contributions of the 
    divergence of a scalar field to the right hand side of the linear equation
    system of the scalar field:
        
    * :math:`A_{ij}=0`
    
    * :math:`A_{ii}=0`
    
    * :math:`b_i=\sum_{f} \dot{m}_f * \phi_f`
    
    * i,j = 1 ... # of cells
    
    .. note::
    
        :math:`\dot{m}_f` is a :class:`pyCFD_fields.fields.SurfaceScalarField`
        value multiplied with 1 if divergence is calculated for the owner cell
        and -1 for the neighbour cell.
        
    .. note::
        :math:`\phi_f` is a :class:`pyCFD_fields.fields.SurfaceScalarField`
        calculated with an iterpolation scheme
        
    .. note::
        boundary face values should be updated before the operator is applied
    
    
    """
    def __init__(self, volume_field, massflux_field, type_):
        r"""
        **Constructor**

        :param volume_field:    volume field to calculate the divergence for.
                                If None than mass divergence is calculated.
        :type volume_field:     pyCFD_fields.fields.VolumeField or None
        :param massflux_field:  surface field with massflux values
        :type massflux_field:   pyCFD_fields.fields.SurfaceScalarField
        :param type_:           type of scheme to calculate the face values 
                                with. Available types are in
                                :class:`pyCFD_fields.calculated_fields.HRSFaceValue` 
                                and :class:`pyCFD_fields.calculated_fields.UpwindFaceValue`.
        :type type_:            string
        """
        mesh_ = massflux_field.father[0]
        pyCFD_operators.generic_operator.GenericScalarOperator.__init__(self, mesh_)
        
        cell_number = len(mesh_.cells)
        if __sparse__:
            self.A         = scipy.sparse.dok_matrix((cell_number, cell_number))
        else:
            self.A         = numpy.zeros((cell_number,cell_number))
        self.b         = numpy.zeros((cell_number,1))
        
        if volume_field is not None:
            if type_ == "UDS":
                face_field = pyCFD_fields.calculated_fields.UpwindFaceValue(volume_field, massflux_field)
            elif type_ == "STOIC":
                face_field = pyCFD_fields.calculated_fields.HRSFaceValue(volume_field, massflux_field, "STOIC")
            elif type_ == "MINMOD":
                face_field = pyCFD_fields.calculated_fields.HRSFaceValue(volume_field, massflux_field, "MINMOD")
            else:
                print type_+" surface interpolation is not supported, use 'UDS', 'MINMOD' or 'STOIC'"
                sys.exit()
            
            for cell_i,cell_ in enumerate(mesh_.cells):
                for face_ in cell_.faces:
                    self.b[cell_i] -= massflux_field.A[face_.id] * face_field.A[face_.id] * face_.get_Sf_sign(cell_)
        else:
            for cell_i,cell_ in enumerate(mesh_.cells):
                for face_ in cell_.faces:
                    self.b[cell_i] -= massflux_field.A[face_.id] * face_.get_Sf_sign(cell_)
            
class DivergenceVec(pyCFD_operators.generic_operator.GenericVectorOperator):
    r"""
    Explicit divergence operator for :class:`pyCFD_fields.fields.VectorField`
    using Picard iteration with known mass fluxes. The operator calls
    :class:`pyCFD_operators.explicit_operators.Divergence` for all the
    components of the field.
    """
    def __init__(self, volume_field, massflux_field, type_):
        r"""
        **Constructor**

        :param volume_field:    volume field to calculate the divergence for
        :type volume_field:     pyCFD_fields.fields.VolumeField
        :param massflux_field:  surface field with massflux values
        :type massflux_field:   pyCFD_fields.fields.SurfaceScalarField
        :param type_:           type of scheme to calculate the face values with. Available types are in pyCFD_fields.calculated_fields.HRSFaceValue .
        :type type_:            string
        """
        mesh_ = volume_field.father[0]
        pyCFD_operators.generic_operator.GenericVectorOperator.__init__(self, mesh_)
        
        # x component
        volume_field_x = volume_field.get_component_as_scalar_field(0)
        divergence_x = Divergence(volume_field_x, massflux_field, type_)
        self.AX = divergence_x.A
        self.bX = divergence_x.b
        
        # y component
        volume_field_y = volume_field.get_component_as_scalar_field(1)
        divergence_y = Divergence(volume_field_y, massflux_field, type_)
        self.AY = divergence_y.A
        self.bY = divergence_y.b
        
        # z component
        volume_field_z = volume_field.get_component_as_scalar_field(2)
        divergence_z = Divergence(volume_field_z, massflux_field, type_)
        self.AZ = divergence_z.A
        self.bZ = divergence_z.b
                
class Laplace(pyCFD_operators.generic_operator.GenericScalarOperator):
    r"""
    Explicit laplace operator for :class:`pyCFD_fields.fields.ScalarField`.
    
    .. math::
        
        \int_{V} \nabla \left( \Gamma \nabla \phi \right) dV &= \oint_{A} \left( \Gamma \nabla \phi \right)_f \cdot \vec{dA} \\
                                                               &= \sum_{f} \left( \Gamma \nabla \phi \right)_f \cdot \vec{A}
    
    The explicit laplace operator returns the explicit contributions of the 
    laplace of a scalar field to the right hand side of the linear equation
    system of the scalar field:
        
    * :math:`A_{ij}=0`
    
    * :math:`A_{ii}=0`
    
    * :math:`b_i=\sum_{f} \left( \Gamma \nabla \phi \right)_f \cdot \vec{A}`
    
    * i,j = 1 ... # of cells
    
    .. note::
        
        Gradient at fixed value boundaries is calculated as:
        
        :math:`\nabla \phi = \frac{\phi_b - \phi_O}{d_{Ob}}`
        
    .. note::
        
        Gradient at fixed gradient boundaries directly applied as:
        
        :math:`\nabla \phi = \nabla \phi_b`
    """
    def __init__(self, volume_field, gamma_ = 1.):
        r"""
        **Constructor**

        :param volume_field:    volume field to calculate the divergence for
        :type volume_field:     pyCFD_fields.fields.VolumeField
        :param gamma_:          default: 1.0, constant conductivity
        :type gamma_:           Float
        """
        mesh_ = volume_field.father[0]
        pyCFD_operators.generic_operator.GenericScalarOperator.__init__(self, mesh_)
        
        cell_number = len(mesh_.cells)
        if __sparse__:
            self.A         = scipy.sparse.dok_matrix((cell_number, cell_number))
        else:
            self.A         = numpy.zeros((cell_number,cell_number))
        self.b         = numpy.zeros((cell_number,1))
        
        for cell_i,cell_ in enumerate(mesh_.cells):
            for face_ in cell_.faces:
                bnd_type = volume_field.patches[face_.bndId].type
                if face_.isBnd == False:
                    if cell_.id == face_.cells[0].id:
                        neighbour_id = face_.cells[1].id
                    else:
                        neighbour_id = face_.cells[0].id
                    # d_on = abs(numpy.linalg.norm(numpy.add(face_.cells[0].C, -face_.cells[1].C)))
                    d_on = abs(cy_operators.cy_linalg_norm(numpy.add(face_.cells[0].C, -face_.cells[1].C)))
                    self.b[cell_i] -= gamma_ * (volume_field.V[neighbour_id] - volume_field.V[cell_.id]) / d_on * face_.A
                elif face_.isBnd and bnd_type == "fixedValue":
                    d_of = abs(cy_operators.cy_linalg_norm(numpy.add(face_.cells[0].C, -face_.C)))
                    self.b[cell_i] -= gamma_ * (volume_field.A[face_.id] - volume_field.V[cell_.id]) / d_of * face_.A
                elif face_.isBnd and bnd_type == "fixedGradient":
                    for patch_face_i,patch_face in enumerate(mesh_.patches[face_.bndId].faces):
                        if patch_face.id == face_.id:
                            self.b[cell_i] -= gamma_ * volume_field.patches[face_.bndId].values[patch_face_i] * face_.A
                else:
                    sys.exit("error in pyCFD_operators.explicit_operators.Laplace, stopping...")
                    
class LaplaceVec(pyCFD_operators.generic_operator.GenericVectorOperator):
    r"""
    Explicit laplace operator for :class:`pyCFD_fields.fields.VectorField`. The
    operator calls :class:`pyCFD_operators.explicit_operators.Laplace` for all
    the components of the field.
    """
    def __init__(self, volume_field, gamma_ = 1.):
        r"""
        **Constructor**

        :param volume_field:    volume field to calculate the divergence for
        :type volume_field:     pyCFD_fields.fields.VolumeField
        :param gamma_:          constant conductivity
        :type gamma_:           Float
        """
        mesh_ = volume_field.father[0]
        pyCFD_operators.generic_operator.GenericVectorOperator.__init__(self, mesh_)
        
        # x component
        volume_field_x = volume_field.get_component_as_scalar_field(0)
        laplace_x = Laplace(volume_field_x, gamma_)
        self.AX = laplace_x.A
        self.bX = laplace_x.b
        
        # y component
        volume_field_y = volume_field.get_component_as_scalar_field(1)
        laplace_y = Laplace(volume_field_y, gamma_)
        self.AY = laplace_y.A
        self.bY = laplace_y.b
        
        # z component
        volume_field_z = volume_field.get_component_as_scalar_field(2)
        laplace_z = Laplace(volume_field_z, gamma_)
        self.AZ = laplace_z.A
        self.bZ = laplace_z.b
                
class Gradient(pyCFD_operators.generic_operator.GenericVectorOperator):
    r"""
    Explicit gradient operator for :class:`pyCFD_fields.fields.VectorField`.
    
    .. math::
        
        \int_{V} \Gamma * \nabla \phi dV = \Gamma * \nabla \phi * \Delta V
                                                               
    The explicit gradient operator returns the explicit contributions of the 
    gradient of a scalar field to the right hand sides of the linear equation
    systems of a vector field:
        
    * :math:`A_{ij}=0`
    
    * :math:`A_{ii}=0`
    
    * :math:`b_{i,x}=\Gamma_i * \left(\nabla \phi \right)_i \cdot \vec{i} * \Delta V_i`
    
    * :math:`b_{i,y}=\Gamma_i * \left(\nabla \phi \right)_i \cdot \vec{j} * \Delta V_i`
    
    * :math:`b_{i,z}=\Gamma_i * \left(\nabla \phi \right)_i \cdot \vec{k} * \Delta V_i`
    
    * i,j = 1 ... # of cells
    
    .. note::
    
        :math:`\nabla \phi` is calculated in
        :class:`pyCFD_fields.calculated_fields.GaussCellGradient`
    """
    def __init__(self, volume_scalar_field, multipl_field=None):
        r"""
        **Constructor**

        :param volume_scalar_field:     volume field to calculate the gradient for
        :type volume_scalar_field:      pyCFD_fields.fields.ScalarField
        :param multipl_field:           default: None, field to multiply the
                                        gradient with (e.g.
                                        :math:`\frac{1}{\rho}`)
        :type multipl_field:            pyCFD_fields.fields.ScalarField
        """
        mesh_ = volume_scalar_field.father[0]
        pyCFD_operators.generic_operator.GenericVectorOperator.__init__(self, mesh_)
        
        gradient_field = pyCFD_fields.calculated_fields.GaussCellGradient(volume_scalar_field)
        
        cell_number = len(mesh_.cells)
        if __sparse__:
            self.AX         = scipy.sparse.dok_matrix((cell_number, cell_number))
            self.AY         = scipy.sparse.dok_matrix((cell_number, cell_number))
            self.AZ         = scipy.sparse.dok_matrix((cell_number, cell_number))
        else:
            self.AX         = numpy.zeros((cell_number,cell_number))
            self.AY         = numpy.zeros((cell_number,cell_number))
            self.AZ         = numpy.zeros((cell_number,cell_number))
        self.bX         = numpy.zeros((cell_number,1))
        self.bY         = numpy.zeros((cell_number,1))
        self.bZ         = numpy.zeros((cell_number,1))
        
        if multipl_field is None:
            for cell_i,cell_ in enumerate(mesh_.cells):
                self.bX[cell_i] -= gradient_field.V[cell_i][0] * cell_.V
                self.bY[cell_i] -= gradient_field.V[cell_i][1] * cell_.V
                self.bZ[cell_i] -= gradient_field.V[cell_i][2] * cell_.V
        else:
            for cell_i,cell_ in enumerate(mesh_.cells):
                self.bX[cell_i] -= gradient_field.V[cell_i][0] * cell_.V * multipl_field.V[cell_i]
                self.bY[cell_i] -= gradient_field.V[cell_i][1] * cell_.V * multipl_field.V[cell_i]
                self.bZ[cell_i] -= gradient_field.V[cell_i][2] * cell_.V * multipl_field.V[cell_i]
