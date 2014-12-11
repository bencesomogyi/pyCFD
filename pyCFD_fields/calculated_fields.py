"""
module for calculated variable fields
"""
__author__ = "Bence Somogyi"
__copyright__ = "Copyright 2014"
__version__ = "0.1"
__maintainer__ = "Bence Somogyi"
__email__ = "bencesomogyi@ivt.tugraz.at"
__status__ = "Prototype"

import numpy
import sys
import pyCFD_fields.fields
import pyCFD_operators.explicit_operators
if sys.platform == 'win32':
    import pyCFD_operators.cython_boost_win32.cy_operators as cy_operators
elif sys.platform == 'linux2':
    import pyCFD_operators.cython_boost_linux2.cy_operators as cy_operators
else:
    sys.exit("unknown platform " + sys.platform + " found in geomTools.py, stopping...")

class LinearFaceValue(pyCFD_fields.fields.Field):
    r"""
    linear interpolation of volume fields resulting in a surface field
    
    **Interpolation:**
    
    Linear interpolation is calculated using weight factors for the 
    neighbour cell values:
    
    .. math:: \phi_f = g_O * \phi_O + g_N * \phi_N
    
    , boundary values are taken from the volue field at the boundary.

    Interpolation weights are calculated in :func:`pyCFD_mesh.face.Face.update_weights()`.
    """
    def __init__(self, volume_field):
        """
        **Constructor**
    
        :param volume_field:    scalar field to calculate face values for
        :type volume_field:     pyCFD_fields.fields.ScalarField
        :return:                a field with face values from linear interpolation
        :rtype:                 pyCFD_fields.fields.Field
        """
        pyCFD_fields.fields.Field.__init__(self)
        field_type = volume_field.type
        mesh_ = volume_field.father[0]
        if (field_type != "scalar") and (field_type != "vector"):
            raise TypeError(
            "unsupported field type '{}' for '{}'"
            ).format(self.__class__, volume_field.type)
        scalar_field = False
        if (field_type == "scalar"):
            scalar_field = True
            
        if scalar_field:
            surface_field = pyCFD_fields.fields.SurfaceScalarField(mesh_, volume_field.name)
        else:
            surface_field = pyCFD_fields.fields.SurfaceVectorField(mesh_, volume_field.name)
        
        self.type = surface_field.type
        self.name = volume_field.name+"_linear_interpolate"
        self.father.append(mesh_)
        self.A = surface_field.A
        
        for face_ind,face_ in enumerate(mesh_.faces):
            # keep boundary value at boundary faces, they should be updated already
            if face_.isBnd:
                if scalar_field: # scalar
                    self.A[face_ind] = volume_field.A[face_ind]
                else: # vector
                    for component in range(2):
                        self.A[face_ind][component] = volume_field.A[face_ind][component]
                continue
            # update internal faces
            owner_id     = mesh_.faces[face_ind].cells[0].id
            neighbour_id = mesh_.faces[face_ind].cells[1].id
            weight_o     = mesh_.faces[face_ind].weights[0]
            weight_n     = mesh_.faces[face_ind].weights[1]
            if scalar_field: # scalar
                value_o      = volume_field.V[owner_id]
                value_n      = volume_field.V[neighbour_id]
                self.A[face_ind] = weight_o * value_o + weight_n * value_n
            else: # vector
                for component in range(2):
                    value_o      = volume_field.V[owner_id][component]
                    value_n      = volume_field.V[neighbour_id][component]
                    self.A[face_ind][component] = weight_o * value_o + weight_n * value_n
                    
class UpwindFaceValue(pyCFD_fields.fields.Field):
    """
    upwind interpolation of volume scalar fields resulting in a surface field
    
    **Interpolation:**
    
    Interpolation is decided based on the value of massflux_field at the faces:
        
    * if the massflux is negative the neighbour cell value is taken
    * otherwise the owner cell value is taken
    """
    def __init__(self, volume_field, massflux_field):
        """
        **Constructor**
        
        :param volume_field:    scalar field to calculate face values for
        :type volume_field:     pyCFD_fields.fields.ScalarField
        :param massflux_field:  surface field with massflux values
        :type massflux_field:   pyCFD_fields.fields.SurfaceScalarField
        :return:                a field including face values from upwind interpolation
        :rtype:                 pyCFD_fields.fields.Field
        """
        mesh_ = volume_field.father[0]
        pyCFD_fields.fields.Field.__init__(self)
        field_type = volume_field.type
        if (field_type != "scalar") and (field_type != "vector"):
            raise TypeError(
            "unsupported field type '{}' for '{}'"
            ).format(self.__class__, volume_field.type)
        scalar_field = False
        if (field_type == "scalar"):
            scalar_field = True
            
        if scalar_field:
            surface_field = pyCFD_fields.fields.SurfaceScalarField(mesh_, volume_field.name)
        else:
            surface_field = pyCFD_fields.fields.SurfaceVectorField(mesh_, volume_field.name)
        
        self.type = surface_field.type
        self.name = volume_field.name
        self.father.append(mesh_)
        self.A = surface_field.A
        
        for face_ind,face_ in enumerate(mesh_.faces):
            # keep boundary value at boundary faces
            if face_.isBnd:
                if scalar_field: # scalar
                    self.A[face_ind] = volume_field.A[face_ind]
                else: # vector
                    for component in range(2):
                        self.A[face_ind][component] = volume_field.A[face_ind][component]
                continue
            # update internal faces
            if massflux_field.A[face_ind] < 0.0:
                cell_ind = face_.cells[1].id
            else:
                cell_ind = face_.cells[0].id
            if scalar_field: # scalar
                self.A[face_ind] = volume_field.V[cell_ind]
            else: # vector
                for component in range(2):
                    self.A[face_ind][component] = volume_field.V[cell_ind][component]
                    
class HRSFaceValue(pyCFD_fields.fields.Field):
    r"""
    High Resolution interpolation of volume scalar fields resulting in 
    surface field
    
    **Interpolation:**
    
    * cell gradient of the volume scalar field is calculated using :class:`GaussCellGradient`
      
    * based on the face massflux central and downwind directions are decided
    .. image:: _images/HRS_U_C_D.png
        :width: 200px
        :align: center
        
    * | upwind value is calculated based on the cell gradient and distance of C and D:
      | :math:`\phi_U = \phi_D - \nabla \phi_P \cdot 2 \vec{d_{CD}}`
      
    * | normalized :math:`\phi` is calculated:
      | :math:`\tilde{\phi_C} = \frac{\phi_C - \phi_D}{\phi_D - \phi_U}`
      
    * | normalized face :math:`\phi` is calculated by according function (for 
        UDS :math:`\tilde{\phi_f} = \tilde{\phi_C}`).
      | Available functions:
          
        * :func:`STOIC`
        
        * :func:`MINMOD`
        
    * | from the normalized face value face value is calculated:
      | :math:`\phi_f = \tilde{\phi_f} \left( \phi_D - \phi_U \right) + \phi_U`

    """
    def __init__(self, volume_field, massflux_field, scheme_):
        """
        **Constructor**
        
        :param volume_field:    scalar field to calculate face values for
        :type volume_field:     pyCFD_fields.fields.ScalarField
        :param massflux_field:  surface field with massflux values
        :type massflux_field:   pyCFD_fields.fields.SurfaceScalarField
        :return:                a field including face values from STOIC scheme
        :rtype:                 pyCFD_fields.fields.Field 
        """
        mesh_ = volume_field.father[0]
        pyCFD_fields.fields.Field.__init__(self)
        field_type = volume_field.type
        if (field_type != "scalar") and (field_type != "vector"):
            raise TypeError(
            "unsupported field type '{}' for '{}'"
            ).format(self.__class__, volume_field.type)
        scalar_field = False
        if (field_type == "scalar"):
            scalar_field = True
            
        if scalar_field:
            surface_field = pyCFD_fields.fields.SurfaceScalarField(mesh_, volume_field.name)
        else:
            surface_field = pyCFD_fields.fields.SurfaceVectorField(mesh_, volume_field.name)
        
        self.type = surface_field.type
        self.name = volume_field.name
        self.father.append(mesh_)
        self.A = surface_field.A
        
        # get gradient field
        grad_phi = GaussCellGradient(volume_field)
        
        for face_ind,face_ in enumerate(mesh_.faces):
            # keep boundary value at boundary faces
            if face_.isBnd:
                if scalar_field: # scalar
                    self.A[face_ind] = volume_field.A[face_ind]
                else: # vector
                    for component in range(2):
                        self.A[face_ind][component] = volume_field.A[face_ind][component]
                continue
            # update internal faces
            if massflux_field.A[face_ind] < 0.0:
                cell_ind_c = face_.cells[1].id
                cell_ind_d = face_.cells[0].id
            else:
                cell_ind_c = face_.cells[0].id
                cell_ind_d = face_.cells[1].id
            # calculate normalized value phi_tilda
            phi_c = volume_field.V[cell_ind_c]
            phi_d = volume_field.V[cell_ind_d]
            r_cd = numpy.add(mesh_.cells[cell_ind_d].C, -mesh_.cells[cell_ind_c].C)
            phi_u = phi_d - numpy.dot(grad_phi.V[cell_ind_c],2.*r_cd)
            phi_f = phi_u
            
            # check if high resolution can be applied
            high_resolution_mode = False
            if (phi_d - phi_u) != 0.:
                high_resolution_mode = True
            
            # apply high resolution scheme
            
            if high_resolution_mode:
                phi_tilda_c = (phi_c - phi_u) / (phi_d - phi_u)
                phi_tilda_f = 0.
                if scheme_ == "STOIC":
                    phi_tilda_f = self.STOIC(phi_tilda_c)
                if scheme_ == "MINMOD":
                    phi_tilda_f = self.MINMOD(phi_tilda_c)
                phi_f = phi_tilda_f * (phi_d - phi_u) + phi_u
            
            if scalar_field: # scalar
                self.A[face_ind] = phi_f
            else: # vector
                for component in range(2):
                    print "In HRSFaceValue a vector version is not implemented yet!"
                    self.A[face_ind][component] = phi_f
                   
    def STOIC(self, phi_tilda_c):
        r"""
        the STOIC high resolution scheme
        
        * | for :math:`0 \leq \tilde{\phi_C} \leq 0.2`
          | :math:`\tilde{\phi_f} = 3*\tilde{\phi_C}`
        
        * | for :math:`0.2 < \tilde{\phi_C} \leq 0.5`
          | :math:`\tilde{\phi_f} = 0.5 + 0.5*\tilde{\phi_C}`
        
        * | for :math:`0.5 < \tilde{\phi_C} \leq 5/6`
          | :math:`\tilde{\phi_f} = 3/8 + 0.75*\tilde{\phi_C}`
        
        * | for :math:`5/6 < \tilde{\phi_C} \leq 1`
          | :math:`\tilde{\phi_f} = 1`
        
        * :math:`\tilde{\phi_f} = \tilde{\phi_C}` elsewhere
        """
        phi_tilda_f = 0.
        if   (phi_tilda_c >= 0.) and (phi_tilda_c<=0.2):
            phi_tilda_f = 3. * phi_tilda_c
        elif (phi_tilda_c > 0.2) and (phi_tilda_c<=0.5):
            phi_tilda_f = 0.5 * phi_tilda_c + 0.5
        elif (phi_tilda_c > 0.5) and (phi_tilda_c<=5./6.):
            phi_tilda_f = 0.75 * phi_tilda_c + 3./8.
        elif (phi_tilda_c > 5./6.) and (phi_tilda_c<=1):
            phi_tilda_f = 1.
        else:
            phi_tilda_f = phi_tilda_c
        return phi_tilda_f
        
    def MINMOD(self, phi_tilda_c):
        r"""
        the MINMOD high resolution scheme
        
        * | for :math:`0 \leq \tilde{\phi_C} \leq 0.5`
          | :math:`\tilde{\phi_f} = 3/2*\tilde{\phi_C}`
        
        * | for :math:`0.5 < \tilde{\phi_C} \leq 1`
          | :math:`\tilde{\phi_f} = 0.5 + 0.5*\tilde{\phi_C}`
        
        * :math:`\tilde{\phi_f} = \tilde{\phi_C}` elsewhere
        """
        phi_tilda_f = 0.
        if   (phi_tilda_c >= 0.) and (phi_tilda_c<=0.5):
            phi_tilda_f = 3./2. * phi_tilda_c
        elif (phi_tilda_c > 0.5) and (phi_tilda_c<=1.):
            phi_tilda_f = 0.5 * phi_tilda_c + 0.5
        else:
            phi_tilda_f = phi_tilda_c
        return phi_tilda_f
                    
class GaussCellGradient(pyCFD_fields.fields.Field):
    r"""
    Calculate cell gradient of a scalar field using the Gauss theorem resulting
    in a vector field:
    
    :math:`\int_{V} \nabla \phi dV = \oint_{A} \phi_f \vec{dA} = \sum_{f} \phi_f \vec{A}`
       
    Iterations are performed to correct non-conjunctionality.
    
    **Iteration:**
    
    1. | Face values are guessed assuming that cells are conjunctional:
       | :math:`\phi_{f'} = g_O * \phi_O + g_N * \phi_N`
       | , the subscripts stand for Owner and Neighbour. Interpolation weights
       | :math:`g_O` and :math:`g_N` are calculated in
       | :func:`pyCFD_mesh.face.Face.update_gradient_weights`
    2. | From guessed face values the volume gradient is calculated:
       | :math:`\nabla \phi_O = \frac{1}{V_O} \sum_{nb} \phi_{f'} \vec{S_f}`
    3. | Update face gradient from cell gradient:
       | :math:`\nabla \phi_f = g_O * {\left( \nabla \phi \right)}_O + g_N {\left( \nabla \phi \right)_N}`
    4. | Update the face value using a correction from the :math:`\phi_{f'}` value
       | :math:`\phi_f = \phi_{f'} + \nabla \phi_f \cdot \left(\vec{r_f}-\vec{r_{f'}}\right)`
    5. | Update :math:`\nabla \phi_O`:
       | :math:`\nabla \phi_O = \frac{1}{V_O} \sum_{nb} \phi_f \vec{S_f}`
    6. Repeat from step 3
    
    **Treatment of boundary conditions:**
    
    At boundaries the value :math:`\phi_f` is known, therefore it is not calculated.
    """
    def __init__(self, volume_field, nonConjIters_ = 0):
        """
        **Constructor**
      
        :param volume_field:    scalar field to calculate the gradient for
        :type volume_field:     :class:`pyCFD_fields.fields.ScalarField`
        :param nonConjIters_:   default: 0, number of non-conjunctional 
                                iteration steps
        :type nonConjIters_:    int
        :return:                a field with cell gradient values
        :rtype:                 :class:`pyCFD_fields.fields.VectorField`
        """
        # define number of iterations
        self.nonConjIters = nonConjIters_
        
        mesh_ = volume_field.father[0]
        pyCFD_fields.fields.Field.__init__(self)
        field_type = volume_field.type
        if (field_type != "scalar"):
            raise TypeError(
            "unsupported field type '{}' for '{}'"
            ).format(self.__class__, volume_field.type)

        self.type = "vector"
        self.name = volume_field.name+"_grad"
        self.father.append(mesh_)
        
        face_phi = LinearFaceValue(volume_field)
        grad_phi = pyCFD_fields.fields.VectorField(mesh_, "grad_phi")
        face_grad_phi = pyCFD_fields.fields.SurfaceVectorField(mesh_, "face_grad_phi")
         
        def update_cell_gradient():
            """
            update cell gradient field as in step 2
            """
            for face_ in mesh_.faces:
                grad_phi.V[face_.cells[0].id] += face_phi.A[face_.id] * face_.Sf / face_.cells[0].V
                if face_.isBnd == False:
                    grad_phi.V[face_.cells[1].id] -= face_phi.A[face_.id] * face_.Sf / face_.cells[1].V
                            
        def update_face_gradient():
            """
            update face gradients as in step 3
            """
            for face_ in mesh_.faces:
                if face_.isBnd == True:
                    continue
                cell_o_id = face_.cells[0].id
                cell_n_id = face_.cells[1].id
                face_grad_phi.A[face_.id] = numpy.add(face_.gradWeights[0]*grad_phi.V[cell_o_id], face_.gradWeights[1]*grad_phi.V[cell_n_id])
                        
        def update_face_value():
            """
            correct face values as in step 4
            """
            for face_ in mesh_.faces:
                if face_.isBnd == True:
                    continue
                # value at fixed gradient boundaries are not corrected!
                face_phi.A[face_.id] += numpy.dot(face_grad_phi.A[face_.id], face_.ffToF)

        for iter_ in range(self.nonConjIters):
            update_cell_gradient()
            update_face_gradient()
            update_face_value()
            iter_ += 1
        
        update_cell_gradient()

        self.V = grad_phi.V
        self.A = face_grad_phi.A
        
class GaussFaceGradient(pyCFD_fields.fields.Field):
    r"""
    Calculate face gradient of a scalar using the iterative loop of the 
    :class:`GaussCellGradient` class but resulting in a surface vector field
    with a stencil reduced to the neighbour cells only.
    
    Reduction of stencil:
        
    :math:`\nabla \phi_f = \overline{\nabla \phi_f} - \left( \overline{\nabla \phi_f} \cdot \vec{e_{ON}}  \right) \vec{e_{ON}} + \frac{\phi_N-\phi_O}{|\vec{d_{ON}}|}\vec{e_{ON}}`
    
    , where O refers to the owner cell, while N to the neighbour cell.
    :math:`\vec{d_{ON}}` is the vector pointing from O's centroid to N's
    centroid and :math:`\vec{e_{ON}}` is the unit vector in the same direction.
    """
    def __init__(self, volume_field, nonConjIters_ = 0):
        """
        **Constructor**
          
        :param volume_field:    scalar field to calculate the gradient for
        :type volume_field:     :class:`pyCFD_fields.fields.ScalarField`
        :return:                a field including gradient face values
        :rtype:                 :class:`pyCFD_fields.fields.SurfaceVectorField`
        """
        # define number of iterations
        self.nonConjIters = nonConjIters_
        
        mesh_ = volume_field.father[0]
        pyCFD_fields.fields.Field.__init__(self)
        field_type = volume_field.type
        if (field_type != "scalar"):
            raise TypeError(
            "unsupported field type '{}' for '{}'"
            ).format(self.__class__, volume_field.type)

        self.type = "vector"
        self.name = volume_field.name+"_grad"
        self.father.append(mesh_)
        
        face_phi = LinearFaceValue(volume_field)
        grad_phi = pyCFD_fields.fields.VectorField(mesh_, "grad_phi")
        face_grad_phi = pyCFD_fields.fields.SurfaceVectorField(mesh_, "face_grad_phi")
            
        def update_cell_gradient():
            """
            update cell gradient field as in step 2
            """
            for cell_ in mesh_.cells:
                # face boundary values should be already updated!
                grad_phi.V[cell_.id] = numpy.array([0., 0., 0.])
                for face_ in cell_.faces:
                    grad_phi.V[cell_.id] += face_phi.A[face_.id] * face_.get_Sf(cell_)
                grad_phi.V[cell_.id] /= cell_.V
                            
        def update_face_gradient():
            """
            update face gradients as in step 3
            """
            for face_ in mesh_.faces:
                if face_.isBnd == True:
                    continue
                cell_o_id = face_.cells[0].id
                cell_n_id = face_.cells[1].id
                face_grad_phi.A[face_.id] = numpy.add(face_.gradWeights[0]*grad_phi.V[cell_o_id], face_.gradWeights[1]*grad_phi.V[cell_n_id])
                        
        def update_face_value():
            """
            correct face values as in step 4
            """
            for face_ in mesh_.faces:
                if face_.isBnd == True:
                    # value at fixed gradient boundaries are not corrected!
                    continue
                face_phi.A[face_.id] += numpy.dot(face_grad_phi.A[face_.id], face_.ffToF)

        for iter_ in range(self.nonConjIters):
            update_cell_gradient()
            update_face_gradient()
            update_face_value()
            iter_ += 1

        if self.nonConjIters > -1:
            for face_i,face_ in enumerate(mesh_.faces):
                # apply fixedGradient boundaries
                if face_.isBnd == True and volume_field.patches[face_.bndId].type == "fixedGradient":
                    patch_ = volume_field.patches[face_.bndId]
                    face_grad_phi.A[face_i] = patch_.values[patch_.find_face_id_in_patch(face_)]
                # calculate boundary gradient at fixedValue boundaries
                elif face_.isBnd == True and volume_field.patches[face_.bndId].type == "fixedValue":
                    patch_ = volume_field.patches[face_.bndId]
                    vec_of = numpy.add(face_.C, -face_.cells[0].C)
                    d_of = cy_operators.cy_linalg_norm(vec_of)
                    # d_of = numpy.linalg.norm(vec_of)
                    unit_vec_of = vec_of / d_of
                    face_grad_phi.A[face_i] = (patch_.values[patch_.find_face_id_in_patch(face_)] - volume_field.V[face_.cells[0].id]) / d_of * unit_vec_of                
                # decrease stencil for face gradients, REFERENCE: Darwish: 9.36,
                else:
                    o_cell = face_.cells[0]
                    n_cell = face_.cells[1]
                    vec_on = numpy.add(n_cell.C, -o_cell.C)
                    d_on = cy_operators.cy_linalg_norm(vec_on)
                    # d_on = numpy.linalg.norm(vec_on)
                    unit_vec_on = vec_on / d_on
                    face_grad_phi.A[face_i] -= numpy.dot(face_grad_phi.A[face_i], unit_vec_on) * unit_vec_on
                    face_grad_phi.A[face_i] += (volume_field.V[n_cell.id] - volume_field.V[o_cell.id]) / d_on * unit_vec_on
                    
        self.A = face_grad_phi.A
        
class MassFlux(pyCFD_fields.fields.SurfaceScalarField):
    r"""
    | Calculate the massflux surface field by:
    | :math:`\dot{m} = \rho_f * \vec{u_f} \cdot \vec{S_f}`
    | , where :math:`\rho_f` and :math:`\vec{u_f}` are the values at the face
    from linear interpolation using the :class:`LinearFaceValue` class.
    """
    def __init__(self, velo_field, rho_field):
        """
        **Constructor**
    
        :param velo_field:  vector field with velocity values
        :type velo_field:   pyCFD_fields.fields.VectorField
        :param rho_field:   scalar field with density values
        :type rho_field:    pyCFD_fields.fields.ScalarField
        :return:            a surface field with the mass flux values
        :rtype:             pyCFD_fields.fields.SurfaceScalarField
        """
        mesh_ = velo_field.father[0]
        field_name = "massFlux"
        
        pyCFD_fields.fields.SurfaceScalarField.__init__(self, mesh_, field_name)
        
        face_velo = LinearFaceValue(velo_field)
        face_rho  = LinearFaceValue(rho_field)
        
        for i_face in range(len(self.A)):
            self.A[i_face] = numpy.dot(face_velo.A[i_face]*face_rho.A[i_face], mesh_.faces[i_face].Sf)
        
class Divergence(pyCFD_fields.fields.ScalarField):
    r"""
    returns a scalar field with the divergence of a field
    
    It uses the explicit operator class
    :class:`pyCFD_operators.explicit_operators.Divergence`
    and feeds the right hand side values to the return field's cell values.
    """
    def __init__(self,  volume_field, massflux_field, type_):
        """
        **Constructor**
        
        :param volume_field:    scalar field to calculate the divergence for
        :type volume_field:     pyCFD_fields.fields.VolumeField
        :param massflux_field:  massflux field which transports the volume_field over the faces
        :type massflux_field:   pyCFD_fields.fields.SurfaceScalarField
        :param type_:           scheme to calculate the face value of transporting velocity
        :type type_:            string
        """
        mesh_ = massflux_field.father[0]
        pyCFD_fields.fields.ScalarField.__init__(self, mesh_, "mass_div")
        temp_div_operator = pyCFD_operators.explicit_operators.Divergence(volume_field, massflux_field, type_)
        self.V = -1. * temp_div_operator.b[:,0]
        
