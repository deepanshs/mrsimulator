
.. _theory:

****************************
How does `mrsimulator` work?
****************************

The line-shape simulation in ``mrsimulator`` is based on the concept of
`Symmetry Pathways in Solid-State NMR` by Grandinetti `et. al.` [#f1]_

Introduction to NMR frequency components
========================================

The nuclear magnetic resonance (NMR) frequency, :math:`\Omega(\Theta, i, j)`,
for the :math:`\left|i\right> \rightarrow \left|j\right>` transition, where
:math:`\left|i\right>` and :math:`\left|j\right>` are the eigenstates of the
stationary-state semi-classical Hamiltonian, can be written as a sum of
frequency components,

.. math::
    :label: eq_1

    \Omega(\Theta, i, j) = \sum_k \Omega_k (\Theta, i, j),

where :math:`\Theta` is the sample's lattice spatial orientation described with
the Euler angles :math:`\Theta = \left(\alpha, \beta, \gamma\right)`, and
:math:`\Omega_k` is the frequency component from the :math:`k^\text{th}`
interaction of the stationary-state semi-classical Hamiltonian.


Each frequency component, :math:`\Omega_k (\Theta, i, j)`, is separated into
three parts,

.. math::
    :label: eq_2

    \Omega_k(\Theta, i, j) = \omega_k ~ \Xi_L^{(k)}(\Theta) ~ \xi_L^{(k)}(i, j),

where :math:`\omega_k` is the size of the :math:`k^\text{th}` frequency
component, and :math:`\Xi_L^{(k)}(\Theta)` and :math:`\xi_L^{(k)}(i, j)` are
the sample's spatial orientation and quantized NMR transition functions
corresponding to the :math:`L^\text{th}` rank spatial and spin irreducible
spherical tensors, respectively.

----

The spatial orientation function, :math:`\Xi_L^{(k)}(\Theta)`, in Eq.
:eq:`eq_2`, is defined in the laboratory frame, where the :math:`z`-axis is the
direction of the external magnetic field. This function is the spatial
contribution to the observed frequency component arising from the
rotation of the :math:`L^\text{th}`-rank irreducible tensor,
:math:`\varrho_{L,n}^{(k)}`, from the principal axis system, to the lab frame
via Wigner rotation which follows,

.. math::
    :label: eq_3

    \Xi_L^{(k)}(\Theta) = \sum_{n_0=-L}^L D^L_{n_0,0}(\Theta_0)
                          \sum_{n_1=-L}^L D^L_{n_1,n_0}(\Theta_1) ~ ... ~
                          \sum_{n_i=-L}^L D^L_{n_i,n}(\Theta_i) ~~
                          \varrho_{L,n}^{(k)}.

.. Here, :math:`\varrho_{L,n}^{(k)}` is defined in the principal axis system of
.. the interaction tensor, here generically denoted with
.. :math:`\boldsymbol{\rho}^{(\lambda)}`, and the subscript
.. :math:`n \in [-L, L]`.
.. The relationship between :math:`\boldsymbol{\rho}^{(\lambda)}` and
.. :math:`\varrho_{L,n}^{(k)}` is described in the next section.

Here, the term :math:`D^L_{n_i,n_j}(\Theta)` is the
`Wigner rotation matrix element <https://en.wikipedia.org/wiki/Wigner_D-matrix>`_,
generically denoted as,

.. math::
    :label: eq_4

    D^L_{n_i,n_j}(\Theta) = e^{-i n_i \alpha} d_{n_i, n_j}^L(\beta) e^{-i n_j \gamma},

where :math:`d_{n_i, n_j}^L(\beta)` is Wigner small :math:`d` element.

----

In the case of the single interaction Hamiltonian, that is, in the absence of
cross-terms, ``mrsimulator`` further defines the product of the size of the
:math:`k^\text{th}` frequency component, :math:`\omega_k`, and the
:math:`L^\text{th}`-rank irreducible tensor components, :math:`\varrho_{L,n}^{(k)}`, in
the principal axis system of the interaction tensor,
:math:`\boldsymbol{\rho}^{(\lambda)}`, as the scaled spatial orientation
tensor (sSOT) components,

.. math::
    :label: eq_5

    \varsigma_{L,n}^{(k)} = \omega_k \varrho_{L,n}^{(k)},

of rank :math:`L`, also defined in the principal axis system of the interaction
tensor, :math:`\boldsymbol{\rho}^{(\lambda)}`.
Using Eqs. :eq:`eq_3` and :eq:`eq_5`, we re-express Eq. :eq:`eq_2` as

.. math::
    :label: eq_6

    \Omega_k(\Theta, i, j) =  \sum_{n_0=-L}^L D^L_{n_0,0}(\Theta_0)
                              \sum_{n_1=-L}^L D^L_{n_1,n_0}(\Theta_1) ~ ... ~
                              \sum_{n_i=-L}^L D^L_{n_i,n}(\Theta_i) ~~
                              \varpi_{L, n}^{(k)}(i,j),

where

.. math::
    :label: eq_7

    \varpi_{L, n}^{(k)}(i,j) = \varsigma_{L,n}^{(k)}~~\xi_L^{(k)}(i, j)

is the frequency tensor components (FT) of rank :math:`L`, defined in the principal
axis system of the interaction tensor and corresponds to the
:math:`\left|i\right> \rightarrow \left|j\right>` spin transition.


.. |quad_description| replace:: The parameter :math:`\omega_q` is defined as
      :math:`\omega_q = \frac{2\piC_q}{2I(2I-1)}`, where :math:`C_q` is the quadrupole
      coupling constant, and :math:`I` is the spin quantum number
      of the quadrupole nucleus. The parameters :math:`\eta_q` and :math:`\omega_0` are the
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
..     - The parameter :math:`\varrho_\text{iso}` is the isotropic nuclear
..       shielding.

..       .. cssclass:: table-bordered table-hover centered
..       .. list-table::
..         :widths: 20 20 60
..         :header-rows: 1

..         * - Order, :math:`M`
..           - Rank, :math:`L`
..           - :math:`\mathbf{\mathcal{R}}_{L,n}`
..         * - 1
..           - 0
..           - :math:`\mathcal{R}_{0,0}^{(\sigma)} = \varrho_\text{iso}`

.. _spatial_orientation_table:

.. |SOF| replace:: :math:`\mathbf{\varsigma}_{L,n}^{(k)}`
.. |L| replace:: :math:`L`
.. |Mth| replace:: :math:`M^\mathrm{th}`


Scaled spatial orientation tensor (sSOT) components in PAS, |SOF|
=================================================================

Single nucleus scaled spatial orientation tensor components
-----------------------------------------------------------

Nuclear shielding interaction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The nuclear shielding tensor, :math:`\boldsymbol{\rho}^{(\sigma)}`, is a second
rank reducible tensor which can be decomposed into a sum of the zeroth-rank
isotropic, first-rank anti-symmetric and second-rank traceless symmetric
irreducible spherical tensors.
In the principal axis system, the zeroth-rank, :math:`\rho_{0,0}^{(\sigma)}`
and the second-rank, :math:`\rho_{2,n}^{(\sigma)}`, irreducible tensors follow,

.. math::
    \begin{array}{c c c c}
    \rho_{0,0}^{(\sigma)} = -\sqrt{3} \sigma_\text{iso}, &
    \rho_{2,0}^{(\sigma)} = \sqrt{\frac{3}{2}} \zeta_\sigma, &
    \rho_{2,\pm1}^{(\sigma)} = 0, &
    \rho_{2,\pm2}^{(\sigma)} = - \frac{1}{2}\eta_\sigma \zeta_\sigma,
    \end{array}

where :math:`\sigma_\text{iso}, \zeta_\sigma`, and :math:`\eta_\sigma` are the
isotropic nuclear shielding, shielding anisotropy, and shielding asymmetry of
the site, respectively. The shielding anisotropy, and asymmetry are defined
using Haeberlen notation.

**First-order perturbation**

The size of the frequency component, :math:`\omega_k`, from the first-order
perturbation expansion of Nuclear shielding Hamiltonian is
:math:`\omega_0=-\gamma B_0`, where :math:`\omega_0` is the Larmor angular
frequency of the nucleus, and :math:`\gamma`, :math:`B_0` are the gyromagnetic
ratio of the nucleus and the macroscopic magnetic flux density of the applied
external magnetic field, respectively. The relation between
:math:`\varrho_{L,n}^{(\sigma)}` and :math:`\rho_{L,n}^{(\sigma)}` follows,

.. math::
    \varrho_{0,0}^{(\sigma)} &= -\frac{1}{\sqrt{3}} \rho_{0,0}^{(\sigma)} \\
    \varrho_{2,n}^{(\sigma)} &=\sqrt{\frac{2}{3}} \rho_{2,n}^{(\sigma)}

.. cssclass:: table-bordered table-striped centered

.. list-table:: A list of scaled spatial orientation tensors in the principal
  axis system of the nuclear shielding tensor, |SOF|, from Eq. :eq:`eq_5` of
  rank L resulting from the Mth order perturbation
  expansion of the Nuclear shielding Hamiltonian is presented.
  :widths: 25 25 50
  :header-rows: 1

  * - Order, :math:`M`
    - Rank, :math:`L`
    - :math:`\varsigma_{L,n}^{(k)} = \omega_k\varrho_{L,n}^{(k)}`

  * - 1
    - 0
    - :math:`\varsigma_{0,0}^{(\sigma)} = -\omega_0\sigma_\text{iso}`

  * - 1
    - 2
    - :math:`\varsigma_{2,0}^{(\sigma)} = -\omega_0 \zeta_\sigma`,

      :math:`\varsigma_{2,\pm1}^{(\sigma)} = 0`,

      :math:`\varsigma_{2,\pm2}^{(\sigma)} = \frac{1}{\sqrt{6}} \omega_0\eta_\sigma \zeta_\sigma`


Electric quadrupole interaction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The electric field gradient (efg) tensor, :math:`\boldsymbol{\rho}^{(q)}`, is
also a second-rank tensor, however, unlike the nuclear shielding tensor, the
efg tensor is always a symmetric second-rank irreducible tensor.
In the principal axis system, this tensor is given as,

.. math::
    \begin{array}{c c c}
    \rho_{2,0}^{(q)} = \sqrt{\frac{3}{2}} \zeta_q, &
    \rho_{2,\pm1}^{(q)} = 0, &
    \rho_{2,\pm2}^{(q)} = - \frac{1}{2}\eta_q \zeta_q,
    \end{array}

where :math:`\zeta_q`, and :math:`\eta_q` are the efg tensor anisotropy, and
asymmetry of the site, respectively. The efg anisotropy, and
asymmetry are defined using Haeberlen convention.

**First-order perturbation**

The size of the frequency component from the first-order perturbation expansion
of Electric quadrupole Hamiltonian is :math:`\omega_k = \omega_q`,
where :math:`\omega_q = \frac{6\pi C_q}{2I(2I-1)}` is the quadrupole splitting
angular frequency. Here, :math:`C_q` is the quadrupole coupling constant, and
:math:`I` is the spin quantum number of the quadrupole nucleus.
The relation between :math:`\varrho_{L,n}^{(q)}` and
:math:`\rho_{L,n}^{(q)}` follows,

.. math::
    \varrho_{2,n}^{(q)} = \frac{1}{3\zeta_q} \rho_{2,n}^{(q)}.

**Second-order perturbation**

The size of the frequency component from the second-order perturbation
expansion of Electric quadrupole Hamiltonian is
:math:`\omega_k = \frac{\omega_q^2}{\omega_0}`, where :math:`\omega_0` is
the Larmor angular frequency of the quadrupole nucleus.
The relation between :math:`\varrho_{L,n}^{(qq)}` and
:math:`\rho_{L,n}^{(q)}` follows,

.. math::
    \varrho_{L,n}^{(qq)} = \frac{1}{9\zeta_q^2} \sum_{m=-2}^2
              \left<L~n~|~2~2~m~n-m\right> \rho_{2,m}^{(q)}~\rho_{2,n-m}^{(q)},

where :math:`\left<L~M~|~l_1~l_2~m_1~m_2\right>` is the Clebsch Gordan
coefficient.

.. cssclass:: table-bordered table-striped centered

.. list-table:: A list of scaled spatial orientation tensors in the principal
                axis system of the efg tensor, |SOF|, from Eq. :eq:`eq_5` of
                rank L resulting from the Mth order perturbation expansion
                of the Electric Quadrupole Hamiltonian is presented.
  :widths: 25 25 50
  :header-rows: 1

  * - Order, :math:`M`
    - Rank, :math:`L`
    - :math:`\varsigma_{L,n}^{(k)} = \omega_k\varrho_{L,n}^{(k)}`

  * - 1
    - 2
    - :math:`\varsigma_{2,0}^{(q)} = \frac{1}{\sqrt{6}} \omega_q`,

      :math:`\varsigma_{2,\pm1}^{(q)} = 0`,

      :math:`\varsigma_{2,\pm2}^{(q)} = -\frac{1}{6} \eta_q \omega_q`

  * - 2
    - 0
    - :math:`\varsigma_{0,0}^{(qq)} = \frac{\omega_q^2}{\omega_0} \frac{1}{6\sqrt{5}} \left(\frac{\eta_q^2}{3} + 1 \right)`

  * - 2
    - 2
    - :math:`\varsigma_{2,0}^{(qq)} = \frac{\omega_q^2}{\omega_0} \frac{\sqrt{2}}{6\sqrt{7}} \left(\frac{\eta_q^2}{3} - 1 \right)`,

      :math:`\varsigma_{2,\pm1}^{(qq)} = 0`,

      :math:`\varsigma_{2,\pm2}^{(qq)} = -\frac{\omega_q^2}{\omega_0} \frac{1}{3\sqrt{21}} \eta_q`

  * - 2
    - 4
    - :math:`\varsigma_{4,0}^{(qq)} = \frac{\omega_q^2}{\omega_0} \frac{1}{\sqrt{70}} \left(\frac{\eta_q^2}{18} + 1 \right)`,

      :math:`\varsigma_{4,\pm1}^{(qq)} = 0`,

      :math:`\varsigma_{4,\pm2}^{(qq)} = -\frac{\omega_q^2}{\omega_0} \frac{1}{6\sqrt{7}} \eta_q`,

      :math:`\varsigma_{4,\pm3}^{(qq)} = 0`,

      :math:`\varsigma_{4,\pm4}^{(qq)} = \frac{\omega_q^2}{\omega_0} \frac{1}{36} \eta_q^2`


.. _spin_transition_theory:

Spin transition functions, :math:`\xi_L^{(k)}(i,j)`
===================================================

The spin transition function is typically
manipulated via the coupling of the nuclear magnetic dipole moment with the
oscillating external magnetic field from the applied radio-frequency pulse.
Considering the strength of the external magnetic rf field is orders of
magnitude larger than the internal spin-couplings, the manipulation of spin
transition functions are described using the orthogonal rotation subgroups.

Single nucleus spin transition functions
----------------------------------------

.. cssclass:: table-bordered table-striped centered

.. list-table:: A list of single nucleus spin transition functions,
                :math:`\xi_L^{(k)}(i,j)`.
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
operators of rank :math:`L`, and :math:`k \in [-L, L]`.
In terms of the tensor product of the Cartesian operators, the above spherical
tensors are expressed as follows,

.. cssclass:: table-bordered table-striped centered

.. list-table::
  :widths: 50 50
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

where :math:`I` is the spin quantum number of the nucleus and
:math:`\hat{\bf{1}}` is the identity operator.

.. cssclass:: table-bordered table-striped centered
.. list-table:: A list of composite single nucleus spin transition functions,
                :math:`\xi_L^{(k)}(i,j)`. Here, `I` is the spin quantum
                number of the nucleus.
  :widths: 50 50
  :header-rows: 1

  * - :math:`\xi_L^{(k)}(i,j)`
    - Value

  * - :math:`\mathbb{c}_0(i,j)`
    - :math:`\frac{4}{\sqrt{125}} \left[I(I+1) - \frac{3}{4}\right] \mathbb{p}(i, j) + \sqrt{\frac{18}{25}} \mathbb{f}(i, j)`

  * - :math:`\mathbb{c}_2(i,j)`
    - :math:`\sqrt{\frac{2}{175}} \left[I(I+1) - \frac{3}{4}\right] \mathbb{p}(i, j) - \frac{6}{\sqrt{35}} \mathbb{f}(i, j)`

  * - :math:`\mathbb{c}_4(i,j)`
    - :math:`-\sqrt{\frac{18}{875}} \left[I(I+1) - \frac{3}{4}\right] \mathbb{p}(i, j) - \frac{17}{\sqrt{175}} \mathbb{f}(i, j)`

.. _frequency_tensor_theory:

Frequency tensor components (FT) in PAS, :math:`\varpi_{L, n}^{(k)}(i,j)`
=========================================================================

Single nucleus frequency tensor components
------------------------------------------

.. cssclass:: table-bordered table-striped centered

.. list-table:: The table presents a list of frequency tensors defined in the principal
  axis system of the respective interaction tensor from Eq. :eq:`eq_7`,
  :math:`\varpi_{L,n}^{(k)}(i,j)`, of rank L resulting from the Mth order perturbation expansion of the interaction Hamiltonians supported in ``mrsimulator``.
  :widths: 20 15 15 50
  :header-rows: 1

  * - Interaction
    - Order, :math:`M`
    - Rank, :math:`L`
    - :math:`\varpi_{L,n}^{(k)}(i,j)`

  * - Nuclear shielding
    - 1
    - 0
    - :math:`\varpi_{0,0}^{(\sigma)}(i,j) = \varsigma_{0,0}^{(\sigma)} ~~ \mathbb{p}(i, j)`

  * - Nuclear shielding
    - 1
    - 2
    - :math:`\varpi_{2,n}^{(\sigma)}(i,j) = \varsigma_{2,n}^{(\sigma)} ~~ \mathbb{p}(i, j)`

  * - Electric Quadrupole
    - 1
    - 2
    - :math:`\varpi_{2,n}^{(q)}(i,j) = \varsigma_{2,n}^{(q)} ~~ \mathbb{d}(i, j)`

  * - Electric Quadrupole
    - 2
    - 0
    - :math:`\varpi_{0,0}^{(qq)}(i,j) = \varsigma_{0,0}^{(qq)} ~~ \mathbb{c}_0(i, j)`

  * - Electric Quadrupole
    - 2
    - 2
    - :math:`\varpi_{2,n}^{(qq)}(i,j) = \varsigma_{2,n}^{(qq)} ~~ \mathbb{c}_2(i, j)`

  * - Electric Quadrupole
    - 2
    - 4
    - :math:`\varpi_{4,n}^{(qq)}(i,j) = \varsigma_{4,n}^{(qq)} ~~ \mathbb{c}_4(i, j)`


.. [#f1] Grandinetti, P. J., Ash, J. T., Trease, N. M. Symmetry pathways in solid-state
    NMR, PNMRS 2011 **59**, *2*, 121-196.
    `DOI: 10.1016/j.pnmrs.2010.11.003 <https://doi.org/10.1016/j.pnmrs.2010.11.003>`_
