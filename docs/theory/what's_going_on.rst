
.. _theory:

============================
How does `mrsimulator` work?
============================

The lineshape simulation in ``mrsimulator`` is build up on the concept of
`Symmetry Pathways in Solid-State NMR <https://www.sciencedirect.com/
science/article/pii/S0079656510001135?via%3Dihub>`_.

NMR frequency components
------------------------

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
    \Omega_k(\Theta, i, j) = \omega_k ~ \Xi_L^{(k)}(\Theta) ~ \xi_L^{(k)}(i, j),
    :label: eq_2

where :math:`\omega_k` is the size of the :math:`k^\text{th}` frequency
component, and :math:`\Xi_L^{(k)}(\Theta)` and :math:`\xi_L^{(k)}(i, j)` are
the sample's spatial orientation and quantized NMR transition functions
corresponding to the :math:`L^\text{th}` rank spatial and spin irreducible
spherical tensors, respectively.

.. _spatial_orientation_theory:

Spatial orientation part
^^^^^^^^^^^^^^^^^^^^^^^^

The spatial orientation function, :math:`\Xi_L^{(k)}(\Theta)`, in Eq.
:eq:`eq_2`, is defined in the laboratory frame where the :math:`z`-axis is the
direction of the external magnetic field.
This function is related to the :math:`L^\text{th}`-rank spatial irreducible
spherical tensors of the interaction in the principal axis system,
:math:`\rho_{L,n_2}^{(k)}`, via Wigner rotation which follows,

.. math::
    \Xi_L^{(k)}(\Theta) = \sum_{n=-L}^L D^L_{n,0}(\Theta) \rho_{L,n}^{(k)},
    :label: eq_3

where :math:`D^L_{n',n}(\Theta)` is the Wigner rotation matrix element,
generically given as,

.. math::
    D^L_{n',n}(\Theta) = e^{-i n' \alpha} d_{n', n}^L(\beta) e^{-i n \gamma}.
    :label: eq_4

Here, :math:`d_{n', n}^L(\beta)` is Wigner small :math:`d` element.

----

In ``mrsimulator``, we further define the product of the size of the
:math:`k^\text{th}` frequency component, :math:`\omega_k`, and the principal
axis system spatial irreducible spherical tensors, :math:`\rho_{L,n_2}^{(k)}`,
of the :math:`k^\text{th}` interaction as

.. math::
    \mathcal{R}_{L,n}^{(k)} = \omega_k~ \rho_{L,n}^{(k)},
    :label: eq_5

and re-express Eq. :eq:`eq_2` as

.. math::
    \Omega_k(\Theta, i, j) = \sum_{n=-L}^L D^L_{n,0}(\Theta)
                               \underbrace{
                                \mathcal{R}_{L,n}^{(k)} ~~ \xi_L^{(k)}(i, j)
                            }_{\text{frequency component}},
    :label: eq_6

where the product, :math:`\mathcal{R}_{L,n}^{(k)} \xi_L^{(k)}(i, j)`, are the
frequency components in the principal axis system corresponding to the
:math:`L^\text{th}` rank spatial irreducible spherical tensors.

.. cssclass:: table-bordered table-hover centered

.. list-table:: A list of supported interactions in ``mrsimulator`` along with
                the corresponding spatial irreducible spherical tensors in the
                principal axis system.
  :widths: 20 10 10 60
  :header-rows: 1

  * - Interaction
    - Rank, :math:`L`
    - .. math::
            \mathbf{\mathcal{R}}_{L,n}
    - Description

  * - First order Nuclear shielding
    - 0
    - .. math::
        \mathcal{R}_{0,0}^{(\sigma)} = \sigma_\text{iso}
    - The parameter :math:`\sigma_\text{iso}` is the isotropic nuclear shielding.

  * - First order Nuclear shielding
    - 2
    - .. math::
        \begin{array}{r l}
          \small
          \mathcal{R}_{2,0}^{(\sigma)} &= \zeta_\sigma, \\
          \mathcal{R}_{2,\pm1}^{(\sigma)} &= 0, \\
          \mathcal{R}_{2,\pm2}^{(\sigma)} &= -\frac{1}{\sqrt{6}} \eta_\sigma \zeta_\sigma
        \end{array}
    - The parameters :math:`\zeta_\sigma` and :math:`\eta_\sigma` are nuclear shielding
      anisotropy and asymmetry parameters defined using Haeberlen convention.

  * - First order Electric Quadrupole
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

  * - Second order Electric Quadrupole
    - 0
    - .. math::
        \mathcal{R}_{0,0}^{(q)} = \frac{\nu_q^2}{\nu_0} \frac{1}{6\sqrt{5}}
            \left(\frac{\eta_q^2}{3} + 1 \right)
    - The parameter :math:`\nu_q` is defined as :math:`\nu_q = \frac{3C_q}{2I(2I-1)}`, where
      :math:`C_q` is the quadrupole coupling constant, and :math:`I` is the spin quantum number
      of the quadrupole nucleus. The parameters :math:`\eta_q` and :math:`\nu_0` are the
      quadrupole asymmetry and Larmor frequency of the nucleus, respectively.

  * - Second order Electric Quadrupole
    - 2
    - .. math::
        \begin{align}
          \mathcal{R}_{2,0}^{(q)} &= \frac{\nu_q^2}{\nu_0} \frac{\sqrt{2}}{6\sqrt{7}}
          \left(\frac{\eta_q^2}{3} - 1 \right), \\
          \mathcal{R}_{2,\pm1}^{(q)} &= 0, \\
          \mathcal{R}_{2,\pm2}^{(q)} &= -\frac{\nu_q^2}{\nu_0} \frac{1}{3\sqrt{21}} \eta_q
        \end{align}
    - The parameter :math:`\nu_q` is defined as :math:`\nu_q = \frac{3C_q}{2I(2I-1)}`, where
      :math:`C_q` is the quadrupole coupling constant, and :math:`I` is the spin quantum number
      of the quadrupole nucleus. The parameters :math:`\eta_q` and :math:`\nu_0` are the
      quadrupole asymmetry and Larmor frequency of the nucleus, respectively.

  * - Second order Electric Quadrupole
    - 4
    - .. math::
        \begin{align}
          \mathcal{R}_{4,0}^{(q)} &= \frac{\nu_q^2}{\nu_0} \frac{1}{\sqrt{70}}
           \left(\frac{\eta_q^2}{18} + 1 \right), \\
          \mathcal{R}_{4,\pm1}^{(q)} &= 0, \\
          \mathcal{R}_{4,\pm2}^{(q)} &= -\frac{\nu_q^2}{\nu_0} \frac{1}{6\sqrt{7}} \eta_q, \\
          \mathcal{R}_{4,\pm3}^{(q)} &= 0, \\
          \mathcal{R}_{4,\pm4}^{(q)} &= \frac{\nu_q^2}{\nu_0} \frac{1}{36} \eta_q^2
        \end{align}
    - The parameter :math:`\nu_q` is defined as :math:`\nu_q = \frac{3C_q}{2I(2I-1)}`, where
      :math:`C_q` is the quadrupole coupling constant, and :math:`I` is the spin quantum number
      of the quadrupole nucleus. The parameters :math:`\eta_q` and :math:`\nu_0` are the
      quadrupole asymmetry and Larmor frequency of the nucleus, respectively.


.. _spin_transition_theory:

Spin transition part
^^^^^^^^^^^^^^^^^^^^

The spin transition function, :math:`\xi_L^{(k)}(i,j)`, is typically
manipulated via the coupling of the nuclear magnetic dipole moment with the
oscillating external magnetic field from the applied radio-frequency pulse.
Considering the strength of the external magnetic rf field is orders of
magnitude larger than the internal spin-couplings, the manipulation of spin
transition functions are described using the orthogonal rotation subgroups
listed here using the lower-case symbols, :math:`\mathbb{s}`,
:math:`\mathbb{p}`, :math:`\mathbb{d}`, and :math:`\mathbb{f}`.

.. cssclass:: table-bordered table-hover centered

.. list-table:: The spin transition function, :math:`\xi_L^{(k)}(i,j)`.
  :widths: 10 15 40 35
  :header-rows: 1

  * - :math:`\xi_L^{(k)}(i,j)`
    - Rank, :math:`L`
    - Value
    - Description

  * - :math:`\mathbb{s}`
    - 0
    - :math:`0`
    - :math:`\left< j | \hat{T}_{00} | j \right> - \left< i | \hat{T}_{00} | i \right>`

  * - :math:`\mathbb{p}`
    - 1
    - :math:`j-i`
    - :math:`\left< j | \hat{T}_{10} | j \right> - \left< i | \hat{T}_{10} | i \right>`

  * - :math:`\mathbb{d}`
    - 2
    - :math:`\sqrt{\frac{3}{2}} \left(j^2 - i^2 \right)`
    - :math:`\left< j | \hat{T}_{20} | j \right> - \left< i | \hat{T}_{20} | i \right>`

  * - :math:`\mathbb{f}`
    - 3
    - :math:`\frac{1}{\sqrt{10}} [5(j^3 - i^3) + (1 - 3I(I+1))(j-i)]`
    - :math:`\left< j | \hat{T}_{30} | j \right> - \left< i | \hat{T}_{30} | i \right>`

.. _irreducible_tensors:

Here, :math:`\hat{T}_{L,k}(\bf{I})` are the irreducible spherical tensor
operators, where :math:`L` is the rank of the tensor and :math:`k \in [-L, L]`.
In terms of the tensor product of the Cartesian operators, the spherical
tensors are expressed as follows.

.. cssclass:: table-bordered table-hover

.. list-table:: We list the spherical tensors up to rank 3 for :math:`k=0`.
                Here, :math:`I` is the spin-quantum number and
                :math:`\hat{\bf{1}}` is the identity operator.
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
