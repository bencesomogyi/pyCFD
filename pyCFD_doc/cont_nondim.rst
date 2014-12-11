.. _nondim:

The non-dimensional equations
=============================

Continuity
----------

The continuity equation for incompressible flow in 3D can be written using
Einstein’s summation as:

.. math::
    :label: conti

    \frac{\partial u_i}{\partial x_i} = 0

As a first step the non dimensional coordinates shall be defined:

.. math::
    :label: nondimdist

    x_i^* := \frac{x_i}{H}

, followed by the non dimensional velocities:

.. math::
    :label: nondimvelo

    u_i^* := \frac{u_i}{U_{\infty}}

Substituting :eq:`nondimdist` and :eq:`nondimvelo` into :eq:`conti` leads to:

.. math::
    :label: nondimvelo1

    \frac{U_{\infty}}{H} \cdot \frac{\partial u_i^*}{\partial x_i^*} = 0

, which can be simplified to:

.. math::
    :label: nondimvelo2

    \frac{\partial u_i^*}{\partial x_i^*} = 0

.. _pcorreq:
	
The continuity equation in vector form:

.. math::
    :label: contivector

    \nabla \vec{v^*} = 0

, where

.. math::
    :label: vectorcomponents

    \vec{v^*} = \left( \begin{array}{c}
    u^* \\
    v^* \\
    w^* \end{array} \right)

Momentum
--------

:math:`x_i` component
^^^^^^^^^^^^^^^^^^^^^

The xi direction momentum equation for incompressible flow in 3D can be written
as:

.. math::
    :label: momdim

    \frac{\partial u_i}{\partial t} + u_j \frac{\partial u_i}{\partial x_j} =
    -\frac{1}{\rho} \frac{\partial p}{\partial x_i}
    + \nu \left( \frac{\partial^2 u_i}{\partial x_j^2} \right)

Non-dimensional variables for spatial coordinates are defined in
:eq:`nondimdist`. New variables shall be defined for non dimensional time and
pressure:

.. math::
    :label: nondimtime

    t^* := \frac{t \cdot U_{\infty}}{H}



.. math::
    :label: nondimpressure

    p^* := \frac{p}{\rho \cdot U_{\infty}^2}

Substituting :eq:`nondimdist`, :eq:`nondimvelo`, :eq:`nondimtime` and 
:eq:`nondimpressure` into :eq:`momdim`:

.. math::
    :label:

    \frac{U_{\infty} \cdot U_{\infty}}{H} \cdot \frac{\partial u_i^*}
    {\partial t^*} + \frac{U_{\infty} \cdot U_{\infty}}{H} \cdot u_j^*
    \frac{\partial u_i^*}{\partial x_j^*} = \\
    -\frac{\rho \cdot U_{\infty} \cdot U_{\infty}}{H} \cdot \frac{1}{\rho}
    \cdot \frac{\partial p^*}{\partial x^*} + \frac{\nu \cdot U_{\infty}}{H^2}
    \cdot \left( \frac{\partial^2 u_i^*}{\partial x_j^2} \right)

Dividing both sides with :math:`\frac{U_{\infty}^2}{H}`, simplifying with ρ in the
pressure term and substituting :math:`Re = \frac{U_{\infty} \cdot H}{\nu}` into
the viscous term, we get the non-dimensional form of the momentum equation in
the :math:`x_i` direction:

.. math::
    :label: nondimpartial

    \frac{\partial u_i^*}{\partial t^*} + u_j^*\frac{\partial u_i^*}
    {\partial x_j^*} = -\frac{\partial p^*}{\partial x^*} + \frac{1}{Re}
    \frac{\partial^2 u_i^*}{\partial x_j^{*2}}

.. _threedvectorform:
	
3D vector form
^^^^^^^^^^^^^^

Based on :eq:`nondimpartial` the momentum equation can be written in vector
form:

.. math::
    :label: nondimvector

    \frac{\partial \vec{v}^*}{\partial t^*} + \vec{v}^* \cdot \nabla \vec{v}^* =
    -\nabla p^* + \frac{1}{Re} \cdot \Delta\left( \vec{v}^* \right)
