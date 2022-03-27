.. _methods_library_documentation:

===============
Methods Library
===============

The syntax for all library methods follows,

.. code-block:: shell

    lib_method = LibraryMethod(
        channels=["29Si"],  # list of isotopes
        magnetic_flux_density=4.7,  # T
        rotor_angle=57.735 * 3.1415 / 180,  # rad
        rotor_frequency=10000,  # Hz
        spectral_dimensions=[
            SpectralDimension(count=512, spectral_width=50000, reference_offset=10),  # dimension-0
            SpectralDimension(count=256, spectral_width=10000, reference_offset=20),  # dimension-1
        ]
    )

where `LibraryMethod` is a placeholder for the library methods and :numref:`table_generic_method`
describes the attributes of the library method object. A list of included methods follows-

* :ref:`bloch_decay`
* :ref:`bloch_decay_CT`
* :ref:`mq_vas`
* :ref:`st_vas`
* :ref:`ssb2d`

.. _bloch_decay:
.. include:: bloch_decay.rst

.. _bloch_decay_CT:
.. include:: bloch_decay_CT.rst

.. _mq_vas:
.. include:: mq_vas.rst

.. _st_vas:
.. include:: st_vas.rst

.. _ssb2d:
.. include:: ssb2d.rst


.. _table_generic_method:
.. cssclass:: table-bordered
.. list-table:: Attribute description for generic library methods.
  :widths: 25 75
  :header-rows: 1

  * - Keywords
    - Description
  * - channels
    - A list of isotope symbols over which the given method applies.
  * - magnetic_flux_density
    - The macroscopic magnetic flux density, in T, of the applied external magnetic field.
  * - rotor_angle
    - The angle between the sample rotation axis and the applied external magnetic field in radians.
  * - rotor_frequency
    - The sample rotation frequency in Hz.
  * - spectral_dimensions
    - A list of spectral dimensions. The coordinates along each spectral dimension is
      described with the keywords, *count* (:math:`N`), *spectral_width*
      (:math:`\nu_\text{sw}`), and *reference_offset* (:math:`\nu_0`). The
      coordinates are given as,

      .. math::
        \left([0, 1, 2, ... N-1] - \frac{T}{2}\right) \frac{\nu_\text{sw}}{N} + \nu_0

      where :math:`T=N` when :math:`N` is even else :math:`T=N-1`.
