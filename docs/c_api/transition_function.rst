.. _transition_function:

Spin Transition Functions
-------------------------

The spin transition function, :math:`\xi_l(i,j)`.

The spin transition function is typically manipulated via the coupling of the
nuclear magnetic dipole moment with the oscillating external magnetic field
from the applied radio-frequency pulse. Considering the strength of the
external magnetic rf field is orders of magnitude larger than the internal
spin-coupling, the manipulation of spin transition functions,
:math:`\xi_l(i,j)`, can also be described using the orthogonal rotation
subgroups described here using the lower-case symbols,
:math:`\mathbb{s}`, :math:`\mathbb{p}`, :math:`\mathbb{d}`, and
:math:`\mathbb{f}`.


:ref:`transition_function_source`

Single nucleus spin transition functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: p(double, double)

.. doxygenfunction:: d(double, double)

.. doxygenfunction:: f(double, double, double)

.. doxygenfunction:: cL(double *, double *, double *, double, double, double)


Two weakly coupled nuclei spin transition functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: dIS(double, double, double, double)
