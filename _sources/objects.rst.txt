

.. _ojects:

-------
Objects
-------

.. _orientation:

Orientation
+++++++++++

  An `Orientation` is a python
  `dict <https://docs.python.org/3/library/stdtypes.html?highlight=dict#dict>`__
  object which represents the three Euler angles. The key-value pairs of this
  object follow,

  .. list-table::
    :widths: 15 25 60
    :header-rows: 1

    * - key
      - value type
      - value description
    * - ``alpha``
      - A `string <https://docs.python.org/3/library/stdtypes.html#str>`__
        containing a physical quantity
      - The :math:`\alpha` Euler angle. For example, '2.1 rad'. The default value is
        '0 rad'.
    * - ``beta``
      - A `string <https://docs.python.org/3/library/stdtypes.html#str>`__
        containing a physical quantity
      - The :math:`\beta` Euler angle. For example, '23.5 deg'. The default value is
        '0 rad'.
    * - ``gamma``
      - A `string <https://docs.python.org/3/library/stdtypes.html#str>`__
        containing a physical quantity
      - The :math:`\gamma` Euler angle. For example, '.5 rad'. The default value is
        '0 rad'.

  *An example of Orientation object.*

.. doctest::
    :skipif: None is None

    >>> {
    ...     "alpha": "0.5 rad",
    ...     "beta": "0.23 rad",
    ...     "gamma": "2.54 rad"
    ... }

.. - The value is a list of three euler angles, [:math:`alpha`, :math:`beta`
..   and :math:`gamma`]. Each angle is given as a string with a physical
..   quantity representing the angle. Tor example, the orientation may be
..   given as ['15 deg', '0.34 rad', '0 rad]. The default value is
..   ['0 rad', '0 rad', '0 rad']




.. _symmetric_tensor:

SymmetricTensor
+++++++++++++++

  A `SymmetricTensor` is a python
  `dict <https://docs.python.org/3/library/stdtypes.html?highlight=dict#dict>`__
  object which represents a traceless second rank symmetric tensor.
  It is represented by the following key-value pairs.

  .. list-table::
    :widths: 15 25 60
    :header-rows: 1

    * - key
      - value type
      - value description
    * - ``anisotropy``
      - A `string <https://docs.python.org/3/library/stdtypes.html#str>`__
        containing a physical quantity
      - The strength of the anisotropy as calculated using Haeberlen
        convention. The value may be provided as a frequency quantity, for
        example, "4.2 kHz", or as a dimensionless frequency ratio, for example,
        "10 ppm" or "10 µHz/Hz". The default value is '0 ppm'
    * - ``asymmetry``
      - A `float <https://docs.python.org/3/library/functions.html#float>`__
      - The asymmetry parameter of the tensor as calculated using
        Haeberlen convention. For example, 0.75. The default value is 0.
    * - ``orientation``
      - An :ref:`orientation` object
      - The Euler angles for rotating the tensor from the principal axis frame
        to the crystal frame.

  *An example of SymmetricTensor object.*

.. doctest::
    :skipif: None is None

    >>> {
    ...     "anisotropy": "10.3 ppm",
    ...     "asymmetry": 0.5,
    ...     "orientation": {
    ...         "alpha": "0.5 rad",
    ...         "beta": "0.23 rad",
    ...         "gamma": "2.54 rad"
    ...     }
    ... }





.. _site:

Site
++++

  A `Site` object is a python
  `dict <https://docs.python.org/3/library/stdtypes.html?highlight=dict#dict>`__
  object which represents a nuclear site with the following key-value pairs,

  .. list-table::
    :widths: 25 25 50
    :header-rows: 1

    * - key
      - value type
      - value description
    * - ``isotope_symbol``
      - A `string <https://docs.python.org/3/library/stdtypes.html#str>`__
      - The NMR active isotope symbol, for example, '13C'.
        This is a required key.
    * - ``isotropic_chemical_shift``
      - A `string <https://docs.python.org/3/library/stdtypes.html#str>`__
        containing a physical quantity
      - The isotropic chemical shift of the isotope. The value may be provided
        as a frequency quantity, "12.6 Hz", or as a dimensionless frequency ratio,
        "1 mHz/Hz", or equivalently, "1000 ppm". The default value is '0 ppm'.
    * - ``shielding_symmetric``
      - A :ref:`symmetric_tensor` object
      - See the description for the :ref:`symmetric_tensor` object.

  *An example of Site object.*

.. doctest::
    :skipif: None is None

    >>> {
    ...     "isotope_symbol": "13C",
    ...     "isotropic_chemical_shift": "15 ppm",
    ...     "shielding_symmetric": {
    ...         "anisotropy": "10.3 ppm",
    ...         "asymmetry": 0.5,
    ...         "orientation": {
    ...             "alpha": "0.5 rad",
    ...             "beta": "0.23 rad",
    ...             "gamma": "2.54 rad"
    ...         }
    ...     }
    ... }




.. _isotopomer:

Isotopomer
++++++++++

  An `Isotopomer` object is a python
  `dict <https://docs.python.org/3/library/stdtypes.html?highlight=dict#dict>`__
  object which represents an isotopomer.
  In `mrsimulator`, each `isotopomer` is treated as a :math:`n`-coupled spin
  system where :math:`n` is the number of sites in the isotopomer.
  It is recommended that if the sites are uncoupled, it be specified as
  individual isotopomers with a single site object, rather than a single
  isotopomer with multiple sites.

  The key-value pairs of the `Isotopomer` object follow,

  .. list-table::
    :widths: 15 15 70
    :header-rows: 1

    * - key
      - value type
      - value description
    * - ``sites``
      - A `list <https://docs.python.org/3/library/stdtypes.html#list>`__
      - A list of :ref:`site` objects.
    * - ``abundance``
      - A `string <https://docs.python.org/3/library/functions.html#float>`__
      - The abundance of the isotopomer. The abundance is given as
        percent, for example, '25.4 %'. This key-value is useful when
        simulating multiple isotopomers. The default value is '100 %'.

  ..  * - ``coulpings``
  ..    - Not yet implemented.


  *An example of Isotopomer object.*

.. doctest::
    :skipif: None is None

    >>> {
    ...     "sites": [
    ...         {
    ...             "isotope_symbol": "13C",
    ...             "isotropic_chemical_shift": "15 ppm",
    ...             "shielding_symmetric": {
    ...                 "anisotropy": "10.3 ppm",
    ...                 "asymmetry": 0.5,
    ...                 "orientation": {
    ...                     "alpha": "0.5 rad",
    ...                     "beta": "0.23 rad",
    ...                     "gamma": "2.54 rad"
    ...                 }
    ...             }
    ...         }
    ...     ],
    ...     "abundance": "15.3 %"
    ... }




.. _direct_dimension:

DirectDimension
+++++++++++++++

  A `DirectDimension` object is a python
  `dict <https://docs.python.org/3/library/stdtypes.html?highlight=dict#dict>`__
  object with the following key-value pairs.

  .. list-table::
    :widths: 25 25 50
    :header-rows: 1

    * - key
      - value type
      - value description
    * - ``nucleus``
      - A `string <https://docs.python.org/3/library/stdtypes.html#str>`__
      - The isotope symbol of the nuclei. The recorded spectrum a histogram of
        frequencies corresponding to this nuclear isotope. An example may
        be '29Si'.
    * - ``magnetic_flux_density``
      - A `string <https://docs.python.org/3/library/stdtypes.html#str>`__
        containing a physical quantity
      - The strength of the external static magnetic field of the spectrometer,
        for example, '14.1 T'. The default value is '9.4 T'.
    * - ``rotor_frequency``
      - A `string <https://docs.python.org/3/library/stdtypes.html#str>`__
        containing a physical quantity
      - The sample spinning frequency, for example, '10 kHz'. The default value
        is '0 Hz'.
    * - ``rotor_angle``
      - A `string <https://docs.python.org/3/library/stdtypes.html#str>`__
        containing a physical quantity
      - The angle between the sample rotation axis and the external magnetic
        field, for example, ‘90 deg’. The default value is ‘54.735 deg’.
    * - ``number_of_points``
      - An `integer <https://docs.python.org/3.3/library/functions.html#int>`__
      - The number of points used in sampling the spectrum, for example, 8192.
        The default value is 1024.
    * - ``spectral_width``
      - A `string <https://docs.python.org/3/library/stdtypes.html#str>`__
        containing a physical quantity
      - The frequency spectral width over which the spectrum is evaluated,
        for example, '500 kHz'. The default value is '100 kHz'.
    * - ``reference_offset``
      - A `string <https://docs.python.org/3/library/stdtypes.html#str>`__
        containing a physical quantity
      - The reference offset of the spectrum, for example, '1.4587 kHz'.
        The default value is '0 Hz'.

.. Note::
    All physical quantities are specified as strings containing a numerical
    value and a unit.
