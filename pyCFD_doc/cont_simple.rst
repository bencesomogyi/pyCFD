.. _simple:

The SIMPLE algorithm
====================

The equations derived in :ref:`nondim` for momentum:

.. math::
    :label:

    \frac{\partial \vec{v}^*}{\partial t^*} + \vec{v}^* \cdot \nabla \vec{v}^* =
    -\nabla p^* + \frac{1}{Re} \cdot \Delta\left( \vec{v}^* \right)

and continuity:

.. math::
    :label: div

    \nabla \vec{v^*} = 0

does not contain an equation for pressure. It should be deducted therefore here.

Writing the momentum equations discretized in space and time, using coefficients of the
linear equation system:

.. math::
    :label: momdiscr

    a_C \cdot \vec{v^*_{pr,C}} + \sum_f \left( a_f \cdot \vec{v^*_{pr,f}} \right)
    = -V_C \cdot \left( \nabla p_{pr}^* \right) + V_C \cdot S_{mom}

, where the subscript *C* indicates the cell center values, *f* the values at
the faces, *a* are the coefficients of the linear equation system and
:math:`S_{mom}` is the source term in the momentum equation. The superscript * is
used for the non dimensional fields. Writing :eq:`momdiscr` for
:math:`\vec{v^*_{pr,C}}`:

.. math::
    :label: momdiscr

    \vec{v^*_{pr,C}} = - H\left( \vec{v^*_{pr,f}} \right)
    -\frac{V_C}{a_C} \cdot \left( \nabla p_{pr}^* \right) + \frac{V_C}{a_C}
    \cdot S_{mom}

, where:

.. math::
    :label: momvelodiscr

    H\left( \vec{v^*_{pr,f}} \right) = \frac{\sum_f \left( a_f \cdot
    \vec{v^*_{pr,f}} \right)}{a_C}

As the solution results in a velocity field which does not divergence free (
:eq:`div` is not fulfilled) this solution is called a predicted velocity field
(or the solution of the momentum equation the predictor step). We search for a
corrective velocity and corrective pressure field

.. math::
    :label: vfinprcorr

    \vec{v_{final}^*} = \vec{v_{pr}^*} + \vec{v_{corr}^*}

.. math::
    :label: pfinprcorr

    p_{final}^* = p_{pr}^* + p_{corr}^*

, which enforces continuity:

.. math::
    :label: newdivfree

    \nabla \vec{v^*_{final}} = 0

Substituting :eq:`momvelodiscr`, :eq:`vfinprcorr` and :eq:`pfinprcorr` into
:eq:`newdivfree`:

.. math::
    :label: pcorrpre

    \nabla\left( \vec{v_{pr}^*} + \vec{v_{corr}^*} \right) = 
    \nabla\left( - H\left( \vec{v^*_{pr,f}} \right)
    -\frac{V_C}{a_C} \cdot \left( \nabla p_{pr}^* \right) + \frac{V_C}{a_C}
    \cdot S_{mom} \right) \\
    + \nabla\left( - H\left( \vec{v^*_{corr,f}} \right)
    -\frac{V_C}{a_C} \cdot \left( \nabla p_{corr}^* \right) + \frac{V_C}{a_C}
    \cdot S_{mom}\right) = 0

The SIMPLE algorithm neglects the corrective terms except for the pressure
correction gradient term. :eq:`pcorrpre` can now be written for
:math:`p_{corr}^*` with only known terms on the right hand side:

.. math::
    :label: pcorreq

    \Delta p_{corr}^* = \frac{a_C}{V_C} \nabla \vec{v^*_{pr,C}}

As :math:`\nabla \vec{v_{pr}^*} = -\nabla \vec{v_{corr}^*}` the velocity
correction becomes:

.. math::
    :label:

    \vec{v_{corr}^*} = - \frac{V_C}{a_C} \nabla p_{corr}^*

Cell velocities are corrected according to :eq:`vfinprcorr`:

.. math::
    :label:

    \vec{v_{final,C}^*} = \vec{v_{pr,C}^*} - \frac{V_C}{a_C} \nabla p_{corr}^*

Correction of the face are also necessery for providing a divergence free
massflux field in the next timestep. The correction is calculated using the
cell gradient of :math:`p_{corr}` interpolated to the faces as:

.. math::
    :label:

    \vec{v_{final,f}^*} = \vec{v_{pr,f}^*} - \left( \frac{V_C}{a_C} \right)_f
    \left( \nabla p_{corr}^* \right)_f







