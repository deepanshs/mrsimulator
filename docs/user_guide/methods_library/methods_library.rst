.. _methods_library_documentation:

===============
Methods Library
===============

For convenience, mrsimulator offers the following pre-built methods -

* :ref:`bloch_decay`
* :ref:`bloch_decay_CT`
* :ref:`mq_vas`
* :ref:`st_vas`
* :ref:`ssb2d`

An example of the syntax that all library methods follows is shown below.

.. code-block:: python

    from mrsimulator.method import SpectralDimension
    from mrsimulator.method.lib import BlochDecaySpectrum

    lib_method = BlochDecaySpectrum(
        channels=["29Si"],  # list of isotopes
        magnetic_flux_density=4.7,  # T
        rotor_angle=57.735 * 3.1415 / 180,  # rad
        rotor_frequency=10000,  # Hz
        spectral_dimensions=[
            SpectralDimension(count=512, spectral_width=5e4, reference_offset=10),
        ],
    )

where `BlochDecaySpectrum` can be replaced with another library method class. Each method has the
*channels* attribute, which is a list of isotopes probed by the method as well as the
*magnetic_flux_density*, *rotor_angle*, and *rotor_frequency* attributes which define the
global experiment parameters.
See :numref:`table_generic_method` for more details.

The method object also has the *spectral_dimensions* attribute, which contains a list of
SpectralDimension objects defining the spectral grid. A 2D method will have two spectral
dimensions in this list, whereas a 1D method will only have one. See
:numref:`table_generic_spectral_dimension` for the attributes of a SpectralDimension object.

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

.. cssclass:: table-bordered table-striped centered

.. _table_generic_method:
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


.. cssclass:: table-bordered table-striped centered

.. _table_generic_spectral_dimension:
.. list-table:: Spectral dimension attributes for use with library methods.
  :widths: 25 75
  :header-rows: 1

  * - Keywords
    - Description

  * - count
    - An integer representing the number of points in the spectral dimension

  * - spectral_width
    - The spectral width of the spectral dimension in Hz.

  * - reference_offset
    - The reference offset of the spectral dimension in Hz.
