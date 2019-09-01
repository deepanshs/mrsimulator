
.. _theory:

****************************
How does `mrsimulator` work?
****************************

The line-shape simulation in ``mrsimulator`` is based on the concept of
`Symmetry Pathways in Solid-State NMR <https://www.sciencedirect.com/
science/article/pii/S0079656510001135?via%3Dihub>`_.

NMR frequency components
========================

The nuclear magnetic resonance (NMR) frequency, :math:`\Omega(\Theta, i, j)`,
for the :math:`i \rightarrow j` transition between the eigenstates of the
stationary-state semi-classical Hamiltonian, can be written as a sum of
frequency components,

.. math::
    \Omega(\Theta, i, j) = \sum_k \Omega_k (\Theta, i, j),
    :label: eq_1

where :math:`\Theta` is the sample's lattice spatial orientation described with
the Euler angles :math:`\Theta = \left(\alpha, \beta, \gamma\right)`, and
:math:`\Omega_k` is the frequency component from the :math:`k^\text{th}`
interaction of the stationary-state semi-classical Hamiltonian.


Each frequency component, :math:`\Omega_k (\Theta, i, j)`, is separated into
three parts,

.. math::
    \Omega_k(\Theta, i, j) = \nu_k ~ \Xi_L^{(k)}(\Theta) ~ \xi_L^{(k)}(i, j),
    :label: eq_2

where :math:`\nu_k` is the size of the :math:`k^\text{th}` frequency
component, and :math:`\Xi_L^{(k)}(\Theta)` and :math:`\xi_L^{(k)}(i, j)` are
the sample's spatial orientation and quantized NMR transition functions
corresponding to the :math:`L^\text{th}` rank spatial and spin irreducible
spherical tensors, respectively.

----

The spatial orientation function, :math:`\Xi_L^{(k)}(\Theta)`, in Eq.
:eq:`eq_2`, is defined in the laboratory frame where the :math:`z`-axis is the
direction of the external magnetic field.
This function is the spatial contribution to the observed frequency component
arising from the :math:`L^\text{th}`-rank spatial irreducible spherical tensors
of the given spin-interaction. It is related to the corresponding spatial
spherical tensors from the principal axis system via Wigner rotation which
follows,

.. math::
    \Xi_L^{(k)}(\Theta) = \sum_{n_0=-L}^L D^L_{n_0,0}(\Theta_0)
                          \sum_{n_1=-L}^L D^L_{n_1,n_0}(\Theta_1) ~ ... ~
                          \sum_{n_i=-L}^L D^L_{n_i,n}(\Theta_i) ~~
                          \rho_{L,n}^{(k)}.
    :label: eq_3

Here, :math:`\rho_{L,n}^{(k)}` is the :math:`L^\text{th}`-rank spatial
irreducible spherical tensor in the principal axis system where subscript,
:math:`n \in [-L, L]`. The term :math:`D^L_{n',n}(\Theta)` is the Wigner
rotation matrix element, generically given as,

.. math::
    D^L_{n',n}(\Theta) = e^{-i n' \alpha} d_{n', n}^L(\beta) e^{-i n \gamma},
    :label: eq_4

where :math:`d_{n', n}^L(\beta)` is Wigner small :math:`d` element.

----

In ``mrsimulator``, we further define the product of the size of the
:math:`k^\text{th}` frequency component, :math:`\nu_k`, and the spatial
irreducible spherical tensor of rank :math:`L` from the principal
axis system, :math:`\rho_{L,n}^{(k)}`, as the spatial orientation function,

.. math::
    \mathcal{R}_{L,n}^{(k)} = \nu_k \rho_{L,n}^{(k)},
    :label: eq_5

of rank :math:`L` defined in the principal axis system.
Using Eqs. :eq:`eq_3` and :eq:`eq_5`, we re-express Eq. :eq:`eq_2` as

.. math::
    \Omega_k(\Theta, i, j) =  \sum_{n_0=-L}^L D^L_{n_0,0}(\Theta_0)
                              \sum_{n_1=-L}^L D^L_{n_1,n_0}(\Theta_1) ~ ... ~
                              \sum_{n_i=-L}^L D^L_{n_i,n}(\Theta_i) ~~
                              \Lambda_{L, n}^{(k)}(i,j),
    :label: eq_6

where

.. math::
    \Lambda_{L, n}^{(k)}(i,j) = \mathcal{R}_{L,n}^{(k)} ~~ \xi_L^{(k)}(i, j)
    :label: eq_7

is the frequency component function, defined in the principal axis system,
corresponding to the :math:`L^\text{th}` rank spatial irreducible spherical
tensors for the :math:`i \rightarrow j` spin transition.


.. |quad_description| replace:: The parameter :math:`\nu_q` is defined as
      :math:`\nu_q = \frac{3C_q}{2I(2I-1)}`, where :math:`C_q` is the quadrupole
      coupling constant, and :math:`I` is the spin quantum number
      of the quadrupole nucleus. The parameters :math:`\eta_q` and :math:`\nu_0` are the
      quadrupole asymmetry and Larmor frequency of the nucleus, respectively.

.. .. cssclass:: table-bordered table-hover centered

.. .. list-table:: A list of :math:`\mathcal{R}_{L,n}^{(k)}` from Eq. :eq:`eq_5`
..                 of rank :math:`L` given in the principal axis system for the
..                 :math:`M^\text{th}` order perturbation expansion of the
..                 interactions supported in ``mrsimulator``.
..   :widths: 20 80
..   :header-rows: 1

..   * - Interaction
..     - Description

..   * - Nuclear shielding
..     - The parameter :math:`\sigma_\text{iso}` is the isotropic nuclear shielding.

..       .. cssclass:: table-bordered table-hover centered
..       .. list-table::
..         :widths: 20 20 60
..         :header-rows: 1

..         * - Order, :math:`M`
..           - Rank, :math:`L`
..           - :math:`\mathbf{\mathcal{R}}_{L,n}`
..         * - 1
..           - 0
..           - :math:`\mathcal{R}_{0,0}^{(\sigma)} = \sigma_\text{iso}`

.. _spatial_orientation_table:

Spatial orientation functions in PAS, :math:`\mathbf{\mathcal{R}}_{L,n}^{(k)}`
------------------------------------------------------------------------------

Single nucleus spin spatial orientation functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered table-hover centered

.. list-table:: A list of spatial orientation functions,
                :math:`\mathcal{R}_{L,n}^{(k)}`, from Eq. :eq:`eq_5` of rank
                :math:`L,` given in the principal axis
                system, for the :math:`M^\text{th}` order perturbation
                expansion of the interactions supported in ``mrsimulator``.
  :widths: 10 8 8 10 64
  :header-rows: 1

  * - Interaction
    - Order, :math:`M`
    - Rank, :math:`L`
    - .. math::
            \mathbf{\mathcal{R}}_{L,n}^{(k)}
    - Description

  * - Nuclear shielding
    - 1
    - 0
    - .. math::
        \mathcal{R}_{0,0}^{(\sigma)} = \nu_0\sigma_\text{iso}
    - The parameter :math:`\sigma_\text{iso}` is the isotropic nuclear
      shielding and :math:`\nu_0` is the Larmor frequency of the nucleus.

  * - Nuclear shielding
    - 1
    - 2
    - .. math::
        \begin{array}{r l}
          \small
          \mathcal{R}_{2,0}^{(\sigma)} &= \nu_0\zeta_\sigma, \\
          \mathcal{R}_{2,\pm1}^{(\sigma)} &= 0, \\
          \mathcal{R}_{2,\pm2}^{(\sigma)} &= -\frac{1}{\sqrt{6}}
                                              \nu_0\eta_\sigma \zeta_\sigma
        \end{array}
    - The parameters :math:`\zeta_\sigma` and :math:`\eta_\sigma` are nuclear
      shielding anisotropy and asymmetry parameters, respectively, defined
      using Haeberlen convention. The parameter :math:`\nu_0` is the Larmor
      frequency of the nucleus.

  * - Electric Quadrupole
    - 1
    - 2
    - .. math::
        \begin{array}{rl}
          \mathcal{R}_{2,0}^{(q)} &= \frac{1}{\sqrt{6}} \nu_q, \\
          \mathcal{R}_{2,\pm1}^{(q)} &= 0, \\
          \mathcal{R}_{2,\pm2}^{(q)} &= -\frac{1}{6} \eta_q \nu_q
        \end{array}
    - The parameter :math:`\nu_q` is defined as :math:`\nu_q = \frac{3C_q}{2I(2I-1)}`, where
      :math:`C_q` is the quadrupole coupling constant, and :math:`I` is the spin quantum number
      of the quadrupole nucleus. The parameter :math:`\eta_q` is the quadrupole asymmetry.

  * - Electric Quadrupole
    - 2
    - 0
    - .. math::
        \mathcal{R}_{0,0}^{(qq)} = \frac{\nu_q^2}{\nu_0} \frac{1}{6\sqrt{5}}
            \left(\frac{\eta_q^2}{3} + 1 \right)
    - |quad_description|

  * - Electric Quadrupole
    - 2
    - 2
    - .. math::
        \begin{align}
          \mathcal{R}_{2,0}^{(qq)} &= \frac{\nu_q^2}{\nu_0} \frac{\sqrt{2}}{6\sqrt{7}}
          \left(\frac{\eta_q^2}{3} - 1 \right), \\
          \mathcal{R}_{2,\pm1}^{(qq)} &= 0, \\
          \mathcal{R}_{2,\pm2}^{(qq)} &= -\frac{\nu_q^2}{\nu_0} \frac{1}{3\sqrt{21}} \eta_q
        \end{align}
    - |quad_description|

  * - Electric Quadrupole
    - 2
    - 4
    - .. math::
        \begin{align}
          \mathcal{R}_{4,0}^{(qq)} &= \frac{\nu_q^2}{\nu_0} \frac{1}{\sqrt{70}}
           \left(\frac{\eta_q^2}{18} + 1 \right), \\
          \mathcal{R}_{4,\pm1}^{(qq)} &= 0, \\
          \mathcal{R}_{4,\pm2}^{(qq)} &= -\frac{\nu_q^2}{\nu_0} \frac{1}{6\sqrt{7}} \eta_q, \\
          \mathcal{R}_{4,\pm3}^{(qq)} &= 0, \\
          \mathcal{R}_{4,\pm4}^{(qq)} &= \frac{\nu_q^2}{\nu_0} \frac{1}{36} \eta_q^2
        \end{align}
    - |quad_description|


.. _spin_transition_theory:

Spin transition functions, :math:`\xi_L^{(k)}(i,j)`
---------------------------------------------------

The spin transition function is typically
manipulated via the coupling of the nuclear magnetic dipole moment with the
oscillating external magnetic field from the applied radio-frequency pulse.
Considering the strength of the external magnetic rf field is orders of
magnitude larger than the internal spin-couplings, the manipulation of spin
transition functions are described using the orthogonal rotation subgroups.

Single nucleus spin transition functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered table-hover centered

.. list-table:: A list of spin transition functions, :math:`\xi_L^{(k)}(i,j)`.
  :widths: 10 12 43 35
  :header-rows: 1

  * - :math:`\xi_L^{(k)}(i,j)`
    - Rank, :math:`L`
    - Value
    - Description

  * - :math:`\mathbb{s}(i,j)`
    - 0
    - :math:`0`
    - :math:`\left< j | \hat{T}_{00} | j \right> - \left< i | \hat{T}_{00} | i \right>`

  * - :math:`\mathbb{p}(i,j)`
    - 1
    - :math:`j-i`
    - :math:`\left< j | \hat{T}_{10} | j \right> - \left< i | \hat{T}_{10} | i \right>`

  * - :math:`\mathbb{d}(i,j)`
    - 2
    - :math:`\sqrt{\frac{3}{2}} \left(j^2 - i^2 \right)`
    - :math:`\left< j | \hat{T}_{20} | j \right> - \left< i | \hat{T}_{20} | i \right>`

  * - :math:`\mathbb{f}(i,j)`
    - 3
    - :math:`\frac{1}{\sqrt{10}} [5(j^3 - i^3) + (1 - 3I(I+1))(j-i)]`
    - :math:`\left< j | \hat{T}_{30} | j \right> - \left< i | \hat{T}_{30} | i \right>`

.. _irreducible_tensors:

Here, :math:`\hat{T}_{L,k}(\bf{I})` are the irreducible spherical tensor
operators, where :math:`L` is the rank of the tensor and :math:`k \in [-L, L]`.
In terms of the tensor product of the Cartesian operators, the spherical
tensors are expressed as follows,

.. cssclass:: table-bordered table-hover

.. list-table::
  :widths: 45 55
  :header-rows: 1

  * - Spherical tensor operator
    - Representation in Cartesian operators
  * - :math:`\hat{T}_{0,0}(\bf{I})`
    - :math:`\hat{1}`
  * - :math:`\hat{T}_{1,0}(\bf{I})`
    - :math:`\hat{I}_z`
  * - :math:`\hat{T}_{2,0}(\bf{I})`
    - :math:`\frac{1}{\sqrt{6}} \left[3\hat{I}^2_z - I(I+1)\hat{1} \right]`
  * - :math:`\hat{T}_{3,0}(\bf{I})`
    - :math:`\frac{1}{\sqrt{10}} \left[5\hat{I}^3_z + \left(1 - 3I(I+1)\right)\hat{I}_z\right]`

where we only list the spherical tensors up to rank 3 for :math:`k=0`.
Here, :math:`I` is the spin-quantum number and :math:`\hat{\bf{1}}` is the
identity operator.

.. cssclass:: table-bordered table-hover centered
.. list-table:: A list of composite spin transition functions,
                :math:`\xi_L^{(k)}(i,j)`. Here, :math:`I` is the spin quantum
                number.
  :widths: 10 90
  :header-rows: 1

  * - :math:`\xi_L^{(k)}(i,j)`
    - Value

  * - :math:`\mathbb{c}_0(i,j)`
    - :math:`\frac{4}{\sqrt{125}} \left[I(I+1) - \frac{3}{4}\right] \mathbb{p}(i, j) + \sqrt{\frac{18}{25}} \mathbb{f}(i, j)`

  * - :math:`\mathbb{c}_2(i,j)`
    - :math:`\sqrt{\frac{2}{175}} \left[I(I+1) - \frac{3}{4}\right] \mathbb{p}(i, j) - \frac{6}{\sqrt{35}} \mathbb{f}(i, j)`

  * - :math:`\mathbb{c}_4(i,j)`
    - :math:`-\sqrt{\frac{18}{875}} \left[I(I+1) - \frac{3}{4}\right] \mathbb{p}(i, j) - \frac{17}{\sqrt{175}} \mathbb{f}(i, j)`

.. _frequency_component_theory:

Frequency component functions in PAS, :math:`\Lambda_{L, n}^{(k)}(i,j)`
-----------------------------------------------------------------------

Single nucleus spin frequency component functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cssclass:: table-bordered table-hover centered

.. list-table:: A list of frequency component functions,
                :math:`\Lambda_{L,n}^{(k)}(i,j)`, from Eq.
                :eq:`eq_7` of rank :math:`L,` given in the principal axis
                system, for the :math:`M^\text{th}` order perturbation
                expansion of the interactions supported in ``mrsimulator``.
  :widths: 20 15 15 50
  :header-rows: 1

  * - Interaction
    - Order, :math:`M`
    - Rank, :math:`L`
    - .. math::
          \Lambda_{L,n}^{(k)}(i,j)

  * - Nuclear shielding
    - 1
    - 0
    - .. math::
        \Lambda_{0,0}^{(\sigma)}(i,j) = \mathcal{R}_{0,0}^{(\sigma)} ~~ \mathbb{p}(i, j)

  * - Nuclear shielding
    - 1
    - 2
    - .. math::
        \Lambda_{2,n}^{(\sigma)}(i,j) = \mathcal{R}_{2,n}^{(\sigma)} ~~ \mathbb{p}(i, j)

  * - Electric Quadrupole
    - 1
    - 2
    - .. math::
        \Lambda_{2,n}^{(q)}(i,j) = \mathcal{R}_{2,n}^{(q)} ~~ \mathbb{d}(i, j)

  * - Electric Quadrupole
    - 2
    - 0
    - .. math::
        \Lambda_{0,0}^{(qq)}(i,j) = \mathcal{R}_{0,0}^{(qq)} ~~ \mathbb{c}_0(i, j)
  * - Electric Quadrupole
    - 2
    - 2
    - .. math::
        \Lambda_{2,n}^{(qq)}(i,j) = \mathcal{R}_{2,n}^{(qq)} ~~ \mathbb{c}_2(i, j)

  * - Electric Quadrupole
    - 2
    - 4
    - .. math::
        \Lambda_{4,n}^{(qq)}(i,j) = \mathcal{R}_{4,n}^{(qq)} ~~ \mathbb{c}_4(i, j)
