"""
generic module for implicit operators
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
import math
import pyCFD_operators.generic_operator
import pyCFD_fields.calculated_fields
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
          
class DdtEuler(pyCFD_operators.generic_operator.GenericScalarOperator):
    r"""
    Euler transient operator for :class:`pyCFD_fields.fields.ScalarField`:
    
    .. math::
        
        \int_{V} \frac{\partial \phi}{\partial t} dV = \frac{\phi^{n+1} - \phi^n}{\Delta t} \Delta V
        
    The operator returns the followng contributions:
    
    * :math:`A_{ij}=0`
    
    * :math:`A_{ii}=\frac{\Delta V_i}{\Delta T}`
    
    * :math:`b_i=\phi^n_i \frac{\Delta V_i}{\Delta T}`
    
    * i,j = 1 ... # of cells
    """
    def __init__(self, equation, dt):
        r"""
        **Constructor**

        :param equation:        field equation of a
                                :class:`pyCFD_fields.fields.ScalarField` to
                                calculate the time derivative for
        :type equation:         :class:`pyCFD_operators.generic_equation.GenericScalarEquation`
        :param dt:              time step
        :type dt:               float
        """
        mesh_ = equation.field.father[0]
        pyCFD_operators.generic_operator.GenericScalarOperator.__init__(self, mesh_)
        
        cell_number = len(mesh_.cells)
        if __sparse__:
            self.A         = scipy.sparse.dok_matrix((cell_number, cell_number))
        else:
            self.A         = numpy.eye(cell_number)
        self.b         = numpy.zeros((cell_number, 1))
        
        for cell_i,cell_ in enumerate(mesh_.cells):
            self.A[cell_i][cell_i] = cell_.V / dt
            self.b[cell_i] = equation.x_old[cell_i] * cell_.V / dt
          
class DdtEulerVec(pyCFD_operators.generic_operator.GenericVectorOperator):
    r"""
    Euler transient operator for :class:`pyCFD_fields.fields.VectorField`.The
    operator calls :class:`pyCFD_operators.implicit_operators.DdtEuler` for all
    the components of the field.
    """
    def __init__(self, equation, dt):
        r"""
        **Constructor**

        :param volume_field:    field equation of a
                                :class:`pyCFD_fields.fields.VectorField` to
                                calculate the time derivative for
        :type volume_field:     :class:`pyCFD_operators.generic_equation.GenericVectorEquation`
        :param dt:              time step
        :type dt:               float
        """
        mesh_ = equation.field.father[0]
        pyCFD_operators.generic_operator.GenericVectorOperator.__init__(self, mesh_)
        
        cell_number = len(mesh_.cells)
        if __sparse__:
            self.A          = scipy.sparse.dok_matrix((cell_number, cell_number))
        else:
            self.A          = numpy.eye(cell_number)
        self.bX         = numpy.zeros((cell_number, 1))
        self.bY         = numpy.zeros((cell_number, 1))
        self.bZ         = numpy.zeros((cell_number, 1))
        
        for cell_i,cell_ in enumerate(mesh_.cells):
            self.A[cell_i, cell_i] = cell_.V / dt
            self.bX[cell_i] = equation.xX_old[cell_i] * cell_.V / dt
            self.bY[cell_i] = equation.xY_old[cell_i] * cell_.V / dt
            self.bZ[cell_i] = equation.xZ_old[cell_i] * cell_.V / dt
            
#class DdtCrankNicholson(pyCFD_operators.generic_operator.GenericScalarOperator):
#    r"""
#    Explicit Second Order Central transient operator
#    
#    
#    """
#    def __init__(self, equation, dt):
#        r"""
#        constructor for Euler transient operator
#
#        :param volume_field:    volume field to calculate the time derivative for
#        :type volume_field:     pyCFD_fields.fields.VolumeField
#        :param dt:              time step
#        :type dt:               float
#        """
#        mesh_ = equation.field.father[0]
#        pyCFD_operators.generic_operator.GenericScalarOperator.__init__(self, mesh_)
#        
#        cell_number = len(mesh_.cells)
#        if __sparse__:
#            self.A         = scipy.sparse.dok_matrix((cell_number, cell_number))
#        else:
#            self.A         = numpy.eye(cell_number)
#        self.b         = numpy.zeros((cell_number, 1))
#        
#        for cell_i,cell_ in enumerate(mesh_.cells):
#            self.A[cell_i, cell_i] = cell_.V / 2. / dt
#            self.b[cell_i] = equation.x_old_old[cell_i] * cell_.V / 2. / dt
            
#class DdtCrankNicholsonVec(pyCFD_operators.generic_operator.GenericVectorOperator):
#    r"""
#    Explicit Second Order Central transient operator for vector fields.
#    """
#    def __init__(self, equation, dt):
#        r"""
#        constructor for Euler transient vector field operator
#
#        :param volume_field:    volume field to calculate the time derivative for
#        :type volume_field:     pyCFD_fields.fields.VolumeField
#        :param dt:              time step
#        :type dt:               float
#        """
#        mesh_ = equation.field.father[0]
#        pyCFD_operators.generic_operator.GenericVectorOperator.__init__(self, mesh_)
#        
#        cell_number = len(mesh_.cells)
#        if __sparse__:
#            self.A          = scipy.sparse.dok_matrix((cell_number, cell_number))
#        else:
#            self.A          = numpy.eye(cell_number)
#        self.bX         = numpy.zeros((cell_number, 1))
#        self.bY         = numpy.zeros((cell_number, 1))
#        self.bZ         = numpy.zeros((cell_number, 1))
#        
#        for cell_i,cell_ in enumerate(mesh_.cells):
#            self.A[cell_i, cell_i] = cell_.V / 2. / dt
#            self.bX[cell_i] = equation.xX_old[cell_i] * cell_.V / 2. / dt
#            self.bY[cell_i] = equation.xY_old[cell_i] * cell_.V / 2. / dt
#            self.bZ[cell_i] = equation.xZ_old[cell_i] * cell_.V / 2. / dt
            
class Divergence(pyCFD_operators.generic_operator.GenericScalarOperator):
    r"""
    Implicit divergence operator for :class:`pyCFD_fields.fields.ScalarField`
    using Picard iteration with known mass fluxes. Divergence is calculated
    using the Gauss theorem:
    
    .. math::
        
        \int_{V} \nabla \left( \rho \vec{v} \phi \right) dV &= \oint_{A} \left( \rho \vec{v} \phi \right) \cdot \vec{A} \\
                                                               &= \sum_{f} \dot{m}_f * \phi_f
                                                               
    The operator creates the coefficient matrix of divergence using the UDS
    scheme for calculating the face value :math:`\phi_f`. High resolution
    schemes are applied as deferred correction :
        
    .. math::
        
        \sum_{f} \dot{m}_f * \phi_f^{HRS,new} = \sum_{f} \dot{m}_f * \phi_f^{UDS,new} + \sum_{f} \dot{m}_f * \left( \phi_f^{HRS,old} - \phi_f^{UDS,old} \right)
        
    The resulting matrix coefficients and right hand side:
    
    * :math:`A_{ij}=-max\left(-\dot{m}_{f,ij},0\right)`
    
    * :math:`A_{ii}=\sum_j max\left(\dot{m}_{f,ij},0\right)`
    
    * :math:`b_i=-\sum_j \dot{m}_{f,ij} \left( \phi_f^{HRS,old} - \phi_f^{UDS,old} \right)`
    
    * i,j = 1 ... # of cells
    
    , where :math:`f` is the face between cell :math:`i` and cell :math:`j`.
    """
    def __init__(self, volume_field, massflux_field, scheme_):
        r"""
        **Constructor**

        :param volume_field:    volume field to calculate the divergence for
        :type volume_field:     :class:`pyCFD_fields.fields.ScalarField`
        :param massflux_field:  surface field with massflux values
        :type massflux_field:   :class:`pyCFD_fields.fields.SurfaceScalarField`
        :param type_:           type of scheme to calculate the face values with. Available types are in
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
            self.A         = numpy.zeros((cell_number, cell_number))
        self.b         = numpy.zeros((cell_number, 1))
        
        for cell_i,cell_ in enumerate(mesh_.cells):
            for face_ in cell_.faces:
                if face_.isBnd:
                    self.b[cell_i] -= massflux_field.A[face_.id] * volume_field.A[face_.id] * face_.get_Sf_sign(cell_)
                    continue             
                self.A[cell_i, cell_i] += max(massflux_field.A[face_.id]*face_.get_Sf_sign(cell_), 0.)
                if cell_i == face_.cells[0].id:
                    self.A[cell_i, face_.cells[1].id] -= max(-massflux_field.A[face_.id]*face_.get_Sf_sign(cell_), 0.)
                else:
                    self.A[cell_i, face_.cells[0].id] -= max(-massflux_field.A[face_.id]*face_.get_Sf_sign(cell_), 0.)
                    
        if scheme_ != "UDS":
            UDS_face_values = pyCFD_fields.calculated_fields.UpwindFaceValue(volume_field, massflux_field)
            HRS_face_values = pyCFD_fields.calculated_fields.HRSFaceValue(volume_field, massflux_field, scheme_)
            for cell_i,cell_ in enumerate(mesh_.cells):
                for face_ in cell_.faces:
                    if face_.isBnd:
                        continue
                    self.b[cell_i] -= massflux_field.A[face_.id] * HRS_face_values.A[face_.id] * face_.get_Sf_sign(cell_)
                    self.b[cell_i] += massflux_field.A[face_.id] * UDS_face_values.A[face_.id] * face_.get_Sf_sign(cell_)
            
class DivergenceVec(pyCFD_operators.generic_operator.GenericVectorOperator):
    r"""
    Implicit divergence operator for :class:`pyCFD_fields.fields.VectorField`
    using Picard iteration with known mass fluxes. Divergence is calculated
    using the Gauss theorem. The operator calls
    :class:`pyCFD_operators.implicit_operators.Divergence` for all the
    components of the field.
    """
    def __init__(self, volume_field, massflux_field, scheme_):
        r"""
        **Constructor**

        :param volume_field:    volume field to calculate the divergence for
        :type volume_field:     :class:`pyCFD_fields.fields.VectorField`
        :param massflux_field:  surface field with massflux values
        :type massflux_field:   :class:`pyCFD_fields.fields.SurfaceScalarField`
        :param type_:           type of scheme to calculate the face values with.
                                Available types are in
                                :class:`pyCFD_fields.calculated_fields.HRSFaceValue` 
                                and :class:`pyCFD_fields.calculated_fields.UpwindFaceValue`.
        :type type_:            string
        """
        mesh_ = massflux_field.father[0]
        pyCFD_operators.generic_operator.GenericVectorOperator.__init__(self, mesh_)
        
        # x component
        volume_field_x = volume_field.get_component_as_scalar_field(0)
        divergence_x = Divergence(volume_field_x, massflux_field, scheme_)
        self.A  = divergence_x.A
        self.bX = divergence_x.b
        
        # y component
        volume_field_y = volume_field.get_component_as_scalar_field(1)
        divergence_y = Divergence(volume_field_y, massflux_field, scheme_)
        self.bY = divergence_y.b
        
        # z component
        volume_field_z = volume_field.get_component_as_scalar_field(2)
        divergence_z = Divergence(volume_field_z, massflux_field, scheme_)
        self.bZ = divergence_z.b
        
class Laplace(pyCFD_operators.generic_operator.GenericScalarOperator):
    r"""
    Implicit laplace operator for :class:`pyCFD_fields.fields.ScalarField`.
    
    .. math::
        
        \int_{V} \nabla \left( \Gamma \nabla \phi \right) dV &= \oint_{A} \left( \Gamma \nabla \phi \right) \cdot \vec{dA} \\
                                                               &= \sum_{f} \left( \Gamma \nabla \phi \right) \cdot \vec{A}
                                                               
    In an arbitrary mesh accuracy of the gradient is calculation is affected
    by mesh non-orthogonality:
        
    .. image:: _images/non_orto.png
        :width: 300px
        :align: center
        
    In the mesh above the vector :math:`\vec{d_{ON}}` or its unit vector
    :math:`\vec{e_{ON}}` is not parallel to the face normal vector
    :math:`\vec{A}`. The face normal vector can be written as a sum of
    :math:`\vec{E}` (parallel to :math:`\vec{e_{ON}}`) and the difference
    :math:`\vec{T}`.
    
    Applying a non-orthogonal correction for the face normal gradients it can
    be rewritten as a sum of orthogonal and non-orthogonal contribution:
        
    .. math::
        
        \nabla \phi \cdot \vec{A} &= \nabla \phi \cdot \vec{E} + \nabla \phi \cdot \vec{T}
                                  &= E \frac{\phi_N - \phi_O}{d_{ON}} + \nabla \phi \cdot \vec{T}
                                  
    The orthogonal part is treated explicitly, while the non-orthogonal part
    with a deferred correction. This treatment results in the following
    coefficient matrix and right hand side:
    
    * :math:`A_{ij}=\Gamma_{f,ij} \frac{E_{f,ij}}{d_{ON,ij}}`
    
    * :math:`A_{ii}=-\sum_j A_{ij}`
    
    * :math:`b_i=-\sum_f \left( \left(\nabla \phi\right)^{old}_f \cdot \vec{T}_f \right)`
    
    * i,j = 1 ... # of cells
    
    , where :math:`f` is the face between cell :math:`i` and cell :math:`j`.
    
    .. note::
        the term :math:`\nabla \left( \phi \right)^{old}_f` is calculated from the
        known :math:`phi` field using
        :class:`pyCFD_fields.calculated_fields.GaussFaceGradient`.
    
    .. note::
        3 non orthogonal correction types are implemented: "minimal correction"
        , "orthogonal correction" and "over relaxed correction" from which the
        "overrelaxed correction" is the one that is tested.
        
    .. note::
        At fixed value boundaries the coefficient matrix and right hand side
        change the following way:
            
        * :math:`A_{ij}=0`
        
        * :math:`b_i-=\Gamma_{f,ij} \frac{E_{f,ij}}{d_{ON,ij}} * \phi_b`
        
    .. note::
        At fixed gradient boundaries the coefficient matrix and right hand side
        change the following way:
            
        * :math:`A_{ij}=0`
        
        * :math:`b_i-=\Gamma_{f,ij} \nabla \phi_b`
    """
    def __init__(self, volume_field, gamma_ = 1., non_ortho_corr = "", calc_matrix = True):
        r"""
        **Constructor**

        :param volume_field:   volume field to calculate the divergence for
        :type volume_field:    :class:`pyCFD_fields.fields.ScalarField`
        :param gamma_:         default: 1.0, constant conductivity
        :type gamma_:          Float
        :param non_ortho_corr: non orthogonal correction type
        :type non_ortho_corr:  String
        :param calc_matrix:    default: Ture, switch to decide if coefficient
                               matrix should be calculated again
        :type calc_matrix:     bool
        """
        mesh_ = volume_field.father[0]
        pyCFD_operators.generic_operator.GenericScalarOperator.__init__(self, mesh_)
        
        cell_number = len(mesh_.cells)
        if __sparse__:
            self.A         = scipy.sparse.dok_matrix((cell_number, cell_number))
        else:
            self.A         = numpy.zeros((cell_number,cell_number))
        self.b         = numpy.zeros((cell_number,1))
        
#        TODO: change to loop for faces
        if non_ortho_corr != "":
            phi_face_grad = pyCFD_fields.calculated_fields.GaussFaceGradient(volume_field, 1)
            
        for cell_i,cell_ in enumerate(mesh_.cells):
            for face_ in cell_.faces:
                bnd_type = volume_field.patches[face_.bndId].type
                if non_ortho_corr == "":
                    if face_.isBnd == False:
                        # d_on = abs(numpy.linalg.norm(numpy.add(face_.cells[1].C, -face_.cells[0].C)))
                        d_on = abs(cy_operators.cy_linalg_norm(numpy.add(face_.cells[1].C, -face_.cells[0].C)))
                        if calc_matrix:
                            self.A[cell_i, cell_i] -= gamma_ / d_on * face_.A
                            if cell_.id == face_.cells[0].id:
                                self.A[cell_i, face_.cells[1].id] = gamma_ / d_on * face_.A
                            else:
                                self.A[cell_i, face_.cells[0].id] = gamma_ / d_on * face_.A
                    elif face_.isBnd and bnd_type == "fixedValue":
                        # d_of = abs(numpy.linalg.norm(numpy.add(face_.C, -cell_.C)))
                        d_of = abs(cy_operators.cy_linalg_norm(numpy.add(face_.C, -cell_.C)))
                        if calc_matrix:
                            self.A[cell_i, cell_i] -= gamma_ / d_of * face_.A
                        self.b[cell_i] -= gamma_ * volume_field.A[face_.id] / d_of * face_.A
                    elif face_.isBnd and bnd_type == "fixedGradient":
                        patch_face_ids = mesh_.patches[face_.bndId].ids
                        patch_face_i = face_.id - patch_face_ids[0]
                        # patch_face_i = patch_face_ids.index(face_.id)
                        self.b[cell_i] -= gamma_ * volume_field.patches[face_.bndId].values[patch_face_i] * face_.A
                    else:
                        sys.exit("error in pyCFD_operators.explicit_operators.Laplace, stopping...")
                else: # non_ortho_corr is defined
                    if face_.isBnd == False:
                        if cell_.id == face_.cells[0].id:
                            vec_cn = numpy.add(face_.cells[1].C, -face_.cells[0].C)
                        else:
                            vec_cn = numpy.add(face_.cells[0].C, -face_.cells[1].C)
                        # d_cn = numpy.linalg.norm(vec_cn)
                        d_cn = cy_operators.cy_linalg_norm(vec_cn)
                        unit_vec_cn = vec_cn / d_cn
                        if   non_ortho_corr == "MINIMUM": # minimum correction approach
                            vec_E = numpy.dot(unit_vec_cn, face_.Sf * face_.get_Sf_sign(cell_)) * unit_vec_cn
                        elif non_ortho_corr == "ORTHOGONAL": # orthogonal correction approach
                            vec_E = face_.A * unit_vec_cn
                        elif non_ortho_corr == "OVERRELAXED": # over-relaxed approach
                            vec_E = math.pow(face_.A, 2.) / numpy.dot(unit_vec_cn, face_.Sf * face_.get_Sf_sign(cell_)) * unit_vec_cn
                        else:
                            print non_ortho_corr + " non orhogonal correction is not implemented, use 'MINIMUM', 'ORTHOGONAL' or 'OVERRELAXED'"
                            sys.exit()
                        vec_T = numpy.add(face_.Sf * face_.get_Sf_sign(cell_), -vec_E)
                        # E = numpy.linalg.norm(vec_E)
                        E = cy_operators.cy_linalg_norm(vec_E)
                        if calc_matrix:
                            self.A[cell_i, cell_i] -= gamma_ / d_cn * E
                            if cell_.id == face_.cells[0].id:
                                self.A[cell_i, face_.cells[1].id] = gamma_ / d_cn * E
                            else:
                                self.A[cell_i, face_.cells[0].id] = gamma_ / d_cn * E
                        self.b[cell_i] -= gamma_ * numpy.dot(phi_face_grad.A[face_.id], vec_T)
                    elif face_.isBnd and bnd_type == "fixedValue":
                        vec_of = numpy.add(face_.C, -cell_.C)
                        # d_of = numpy.linalg.norm(vec_of)
                        d_of = cy_operators.cy_linalg_norm(vec_of)
                        unit_vec_of = vec_of / d_of
                        if   non_ortho_corr == "MINIMUM": # minimum correction approach
                            vec_E = numpy.dot(unit_vec_of, face_.Sf * face_.get_Sf_sign(cell_)) * unit_vec_of
                        elif non_ortho_corr == "ORTHOGONAL": # orthogonal correction approach
                            vec_E = face_.A * unit_vec_of
                        elif non_ortho_corr == "OVERRELAXED": # over-relaxed approach
                            vec_E = math.pow(face_.A, 2.) / numpy.dot(unit_vec_of, face_.Sf * face_.get_Sf_sign(cell_)) * unit_vec_of
                        else:
                            print non_ortho_corr + " non orhogonal correction is not implemented, use 'MINIMUM', 'ORTHOGONAL' or 'OVERRELAXED'"
                            sys.exit()
                        vec_T = numpy.add(face_.Sf * face_.get_Sf_sign(cell_), -vec_E)
                        # E = numpy.linalg.norm(vec_E)
                        E = cy_operators.cy_linalg_norm(vec_E)
                        if calc_matrix:
                            self.A[cell_i, cell_i] -= gamma_ / d_of * E
                        self.b[cell_i] -= gamma_ * volume_field.A[face_.id] / d_of * E
                        self.b[cell_i] -= gamma_ * numpy.dot(phi_face_grad.A[face_.id], vec_T)
                    elif face_.isBnd and bnd_type == "fixedGradient":
                        for patch_face_i,patch_face in enumerate(mesh_.patches[face_.bndId].faces):
                            if patch_face.id == face_.id:
                                self.b[cell_i] -= gamma_ * volume_field.patches[face_.bndId].values[patch_face_i] * face_.A
                    else:
                        sys.exit("error in pyCFD_operators.explicit_operators.Laplace, stopping...")
        
class LaplaceVec(pyCFD_operators.generic_operator.GenericVectorOperator):
    r"""
    Implicit laplace operator for :class:`pyCFD_fields.fields.VectorField`.
    The operator calls :class:`pyCFD_operators.implicit_operators.Laplace` for
    all the components of the field.
    """
    def __init__(self, volume_field, gamma_ = 1., non_ortho_corr = "", calc_matrix = True):
        r"""
        constructor for Divergence operator

        :param volume_field:                volume field to calculate the divergence for
        :type volume_field:                 pyCFD_fields.fields.VolumeField
        :param gamma_:                      constant conductivity
        :type gamma_:                       Float
        :param gamnon_ortho_corr:           non orthogonal correction type
        :type gamma_:                       String
        """
        mesh_ = volume_field.father[0]
        pyCFD_operators.generic_operator.GenericVectorOperator.__init__(self, mesh_)
        
        # x component
        volume_field_x = volume_field.get_component_as_scalar_field(0)
        laplace_x = Laplace(volume_field_x, gamma_, non_ortho_corr, calc_matrix)
        self.A  = laplace_x.A
        self.bX = laplace_x.b
        
        # y component
        volume_field_y = volume_field.get_component_as_scalar_field(1)
        laplace_y = Laplace(volume_field_y, gamma_, non_ortho_corr, calc_matrix)
        self.bY = laplace_y.b
        
        # z component
        volume_field_z = volume_field.get_component_as_scalar_field(2)
        laplace_z = Laplace(volume_field_z, gamma_, non_ortho_corr, calc_matrix)
        self.bZ = laplace_z.b