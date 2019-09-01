

===============
Spatial Tensors
===============

The NMR frequency can be written as a sum of frequency components where each
component :math:`\nu_k(\Theta, i, j)` is given as a product of the spatial and
spin functions which are determined from the spatial and spin part of the
Hamiltonian, given as,

.. math::
    \nu_k(\Theta, i, j) = \Xi_k(\Theta) \xi_k(i, j),

where :math:`\Xi_k(\Theta)` is the spatial and :math:`\xi_k(i, j)` is the
spin part of the :math:`k^\text{th}` interaction Hamiltonian. Here,
:math:`\Theta` is the
lattice spacial orientation given by three Euler angles, :math:`\alpha`,
:math:`\beta`, and :math:`\gamma`, and :math:`i \rightarrow j` is the
transition between the eigenstates of the stationary-state semi-classical
Hamiltonain.

The spatial part of the spin Hamiltonian, :math:`\Xi(\Theta)`, may be expressed
as irreducible tensor of rank :math:`l`. In the principal axis system, PAS,
this tensor may be expressed in terms of the principal values of the
interaction tensor. In the following, we list the spatial irreducible
tensor(s) in the principal axis system, :math:`\rho_{l,m}`, of rank :math:`l`,
resulting from the interaction Hamiltonian.


First order Nuclear shielding
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: spatial_tensors_from_1st_order_nuclear_shielding_Hamiltonian(double
    *restrict, complex128 *restrict, const double, const double, const double)

First order Electric Quadrupole
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: spatial_tensors_from_1st_order_electric_quadrupole_Hamiltonian(complex128
    *restrict, const double, const double, const double)

Second order Electric Quadrupole
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: spatial_tensors_from_2nd_order_electric_quadrupole_Hamiltonian(double
    *restrict, complex128 *restrict, complex128 *restrict, const double, const
    double, const double, const double, const int)
