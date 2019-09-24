

.. _spectroscopic_dimension:

Understanding SpectroscopicDimension
------------------------------------

We define a `SpectroscopicDimension` as an object with attributes representing
a spectroscopic dimension of an NMR experiment.

.. cssclass:: table-bordered table-hover

.. list-table:: Attributes of a Dimension object.
  :widths: 25 25 50
  :header-rows: 1

  * - Attributes
    - Type
    - Description

  * - ``isotope``
    - String
    - An `optional` isotope string given as the atomic number followed by
      the isotope symbol, for example, ``13C``, ``29Si``, ``27Al``, and so on.

  * - ``magnetic_flux_density``
    - ScalarQuantity
    - The magnetic flux density of the external static magnetic field of the
      spectrometer, for example, ``14.1 T``. The default value is ``9.4 T``.

  * - ``rotor_frequency``
    - ScalarQuantity
    - The sample spinning frequency, for example, ``10 kHz``. The default value
      is ``0 Hz``.

  * - ``rotor_angle``
    - ScalarQuantity
    - The angle between the sample rotation axis and the external magnetic
      field, for example, ``1 rad``. The default value is ``54.735 deg``.

  * - ``number_of_points``
    - Integer
    - The number of points sampled along the spectroscopic dimension, for example, 8192.
      The default value is 1024.

  * - ``spectral_width``
    - ScalarQuantity
    - The frequency spectral width along the spectroscopic dimension over which
      the spectrum is evaluated, for example, ``500 kHz``. The default value
      is ``100 kHz``.

  * - ``reference_offset``
    - ScalarQuantity
    - The reference offset along the spectroscopic dimension, for example,
      ``1.4587 kHz``. The default value is ``0 Hz``.

**Example**

Using python `dict <https://docs.python.org/3/library/stdtypes.html?highlight=dict#dict>`_
object, the `Dimension` object may be represented as,

.. code-block:: py

    spectroscopic_dimension_object = dict(
        "isotope": "13C",
        "magnetic_flux_density": "9.4 T",
        "rotor_frequency": "5 kHz",
        "rotor_angle": "54.735 deg",
        "number_of_points": 8192,
        "spectral_width": "100 kHz",
        "reference_offset": "0 Hz",
    )
