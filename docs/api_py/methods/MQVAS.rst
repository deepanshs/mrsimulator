.. _mqvas_ref:

Multi-quantum variable-angle spinning
-------------------------------------

The following classes are used in simulating a multi-quantum variable-angle spinning spectrum.
The resulting spectrum is sheared and scaled, such that the frequencies along the indirect
dimension is given as

.. math::
    \langle \Omega\rangle_\text{MQ-VAS} =  \frac{1}{1+\kappa}\Omega_{m, -m} +
            \frac{\kappa}{1+\kappa}\Omega_{1/2, -1/2}.

Here, :math:`\langle \Omega\rangle_\text{MQ-VAS}` is the average frequency along the indirect
dimension, :math:`\Omega_{m, -m}` and :math:`\Omega_{1/2, -1/2}` are the frequency
contributions from the :math:`|m\rangle \rightarrow |-m\rangle` symmetric multi-quantum
transition and the central transition, respectively, and :math:`\kappa` is the shear factor.
The values for the shear factor for various satellite transitions are listed below.

.. cssclass:: table-bordered

.. list-table:: The shear factor, :math:`\kappa`, for various symmetric multi-quantum transitions
        in an MQ-VAS method.
    :widths: 30 40 30
    :header-rows: 1

    * - Spin
      - Symmetric multi-quantum transition
      - :math:`\kappa`

    * - 3/2
      - :math:`\left(\frac{3}{2} \rightarrow -\frac{3}{2}\right)`
      - 21/27

    * - 5/2
      - :math:`\left(-\frac{3}{2} \rightarrow \frac{3}{2}\right)`
      - 114/72

    * - 5/2
      - :math:`\left(\frac{5}{2} \rightarrow -\frac{5}{2}\right)`
      - 150/72

    * - 7/2
      - :math:`\left(-\frac{3}{2} \rightarrow \frac{3}{2}\right)`
      - 303/135

    * - 7/2
      - :math:`\left(-\frac{5}{2} \rightarrow \frac{5}{2}\right)`
      - 165/135

    * - 9/2
      - :math:`\left(-\frac{3}{2} \rightarrow \frac{3}{2}\right)`
      - 546/216

    * - 9/2
      - :math:`\left(-\frac{5}{2} \rightarrow \frac{5}{2}\right)`
      - 570/216



.. cssclass:: table-bordered

.. list-table::  The parameters for all MQVAS methods classes.
  :widths: 25 75
  :header-rows: 1

  * - Keywords
    - Description
  * - channels
    - A list of isotope symbols over which the given method applies.
  * - magnetic_flux_density
    - An `optional` float containing the macroscopic magnetic flux density, :math:`H_0`,
      of the applied external magnetic field in units of T. The default value is 9.4.
  * - rotor_angle
    - An `optional` float containing the angle between the sample rotation axis and the
      applied external magnetic field, :math:`\theta`, in units of rads. The default value
      is 0.9553166, i.e. the magic angle.
  * - spectral_dimensions
    - A list of spectral dimensions given as python dictionary object. Each dictionary object
      contains keywords that describe the coordinates along a spectral dimension. The keywords
      along with its definition are:

      - count:
          An optional integer with the number of points, :math:`N`, along the
          dimension. The default value is 1024.
      - spectral_width:
          An `optional` float with the spectral width, :math:`\nu_\text{sw}`, along the
          dimension in units of Hz. The default is 25 kHz.
      - reference_offset:
          An `optional` float with the reference offset, :math:`\nu_0` along the
          dimension in units of Hz. The default value is 0 Hz.
      - origin_offset:
          An `optional` float with the origin offset (Larmor frequency) along the
          dimension in units of Hz. The default value is None.

      The coordinates along each spectral dimension are given as,

      .. math::
        \left([0, 1, 2, ... N-1] - \frac{T}{2}\right) \frac{\nu_\text{sw}}{N} + \nu_0

      where :math:`T=N` when :math:`N` is even else :math:`T=N-1`.

.. note::
    The `rotor_frequency` parameter is fixed for this method. The method produces an
    infinite spinning speed spectrum.


.. currentmodule:: mrsimulator.methods

.. _threeQ_vas_ref:
Triple-quantum variable-angle spinning method
'''''''''''''''''''''''''''''''''''''''''''''

.. autoclass:: ThreeQ_VAS

.. _fiveQ_vas_ref:
Five-quantum variable-angle spinning method
'''''''''''''''''''''''''''''''''''''''''''

.. autoclass:: FiveQ_VAS
