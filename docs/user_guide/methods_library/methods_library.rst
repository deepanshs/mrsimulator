.. _methods_library_documentation:

===============
Methods Library
===============

For convenience, **MRSimulator** offers the following pre-built methods -

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
``channels`` attribute, which is a list of isotopes probed by the method as well as the
``magnetic_flux_density``, ``rotor_angle``, and ``rotor_frequency`` attributes which define the
global experiment parameters.
See :numref:`table_generic_method` for more details.

The method object also has the ``spectral_dimensions`` attribute, which contains a list of
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


Attribute Summaries
-------------------

.. cssclass:: table-bordered table-striped centered
.. _table_generic_method:
.. list-table:: Attribute description for generic library methods.
  :widths: 20 15 65
  :header-rows: 1

  * - Attribute Name
    - Type
    - Description

  * - channels
    - ``list``
    - A *required* list of isotopes given as strings over which the given method applies.
      For example, ``["1H"]``.

  * - magnetic_flux_density
    - ``float``
    - An *optional* float describing the macroscopic magnetic flux density of the applied
      external magnetic field in tesla. For example, ``18.8`` tesla. The default value is
      ``9.4`` tesla.

  * - rotor_frequency
    - ``float``
    - An *optional* float describing the sample rotation frequency in Hz. For example, ``2000`` Hz.
      The default value is ``0`` Hz.

  * - rotor_angle
    - ``float``
    - An *optional* float describing the angle between the sample rotation axis and the external
      magnetic field in radians. The default value is the magic angle,
      ``54.735 * 3.14159 / 180 = 0.955305`` radians.

  * - spectral_dimensions
    - ``list``
    - A list of :ref:`spectral_dim_api` objects describing the spectral dimensions for the method.

  * - simulation
    - CSDM object
    - A CSDM object representing the spectrum simulated by the method. By default, the value is
      ``None``. A value is assigned to this attribute when you run the
      simulation using the :py:meth:`~mrsimulator.Simulator.run` method.

  * - experiment
    - CSDM object
    - An *optional* CSDM object holding an experimental measurement of the method. The default
      value is ``None``


.. cssclass:: table-bordered table-striped centered
.. _table_generic_spectral_dimension:
.. list-table:: Spectral dimension attributes for use with library methods.
  :widths: 20 15 65
  :header-rows: 1

  * - Attribute Name
    - Type
    - Description

  * - count
    - ``int``
    - An *optional* integer representing the number of points, :math:`N`, along the spectroscopic
      dimension. For example, ``4096``. The default value is ``1024``.

  * - spectral_width
    - ``float``
    - An *optional* float representing the width, :math:`\Delta x`, of the spectroscopic dimension
      in Hz. For example, ``10e3`` for 10 kHz. The default value is ``25000`` Hz.

  * - reference_offset
    - ``float``
    - An *optional* float representing the reference offset, :math:`x_0`, of the spectroscopic
      dimension in Hz. For example, ``-8000`` Hz. The default value is ``0``.

  * - origin_offset
    - ``float``
    - An *optional* float representing the origin offset, or Larmor frequency, along the
      spectroscopic dimension in units of Hz. The default value is ``None`` and the origin offset
      is set to the Larmor frequency of isotope from the :attr:`~mrsimulator.Method.channels`
      attribute of the method containing the spectral dimension.
