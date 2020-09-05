Multi-quantum magic angle spinning method
-----------------------------------------

Simulate a multiple-quantum magic-angle spinning spectrum using average frequency. The
resulting spectrum is sheared such that the correlating dimensions are the isotropic
dimension and the MAS dimension, respectively, where the isotropic dimension is given as

.. math::
    \langle \Omega_\text{iso}\rangle_\text{MQ-MAS} =
            \frac{1}{1+\kappa}\Omega_\text{iso}(m, -m) +
            \frac{\kappa}{1+\kappa}\Omega_\text{iso}\left(\frac{1}{2}, -\frac{1}{2}\right),

where :math:`\langle \Omega_\text{iso}\rangle_\text{MQ-MAS}` is the average isotropic
frequency along the isotropic dimension, :math:`\kappa` is the shear factor,
:math:`\Omega_\text{iso}(m, -m)` is the isotropic contribution from a
:math:`|m\rangle \rightarrow |-m\rangle` symmetric multi-quantum transition, and
:math:`\Omega_\text{iso}\left(\frac{1}{2}, -\frac{1}{2}\right)` is the isotropic contributions
from the central transition. The values for the shear factor for various symmetric transitions
are listed below.

.. list-table::
    :widths: 20 40 40
    :header-rows: 1

    * - Spin
      - :math:`\kappa` for transition :math:`\left(\frac{3}{2} \rightarrow -\frac{3}{2}\right)`
      - :math:`\kappa` for transition :math:`\left(\frac{5}{2} \rightarrow -\frac{5}{2}\right)`

    * - 3/2
      - 21/27
      - -

    * - 5/2
      - 114/72
      - 150/72

    * - 7/2
      - 303/135
      - 165/135

    * - 9/2
      - 546/216
      - 570/216

.. currentmodule:: mrsimulator.methods

Triple-quantum magic angle spinning method
''''''''''''''''''''''''''''''''''''''''''

.. currentmodule:: mrsimulator.methods

.. autoclass:: ThreeQ_VAS
.. autoclass:: FiveQ_VAS
