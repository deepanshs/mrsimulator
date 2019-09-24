

.. _dictionary_objects:

*************************
Understanding Isotopomers
*************************

Isotopomers are collection of isolated spin-systems where each isolated
spin-system is a collection of sites and couplings. We refer an isolated
spin-system as an isotopomer. In our model, each isotopomer is composed
of various objects. In the following, we list and describe these objects.

.. A Unified Modeling Language (UML) class-diagram for `Isotopomers` object is
.. shown below.

.. .. figure:: ./_static/class-diagram.png

.. _orientation:

Orientation
-----------

We define `Orientation` as an object with attributes representing the Euler
angles, :math:`\alpha, \beta, \gamma`, involved in the rotation of the
interaction tensor from the principal axis system (PAS) to a reference frame.


.. cssclass:: table-bordered table-hover
.. list-table:: Attributes of an Orientation object.
  :widths: 15 20 65
  :header-rows: 1

  * - Attribute name
    - Type
    - Description

  * - ``alpha``
    - ScalarQuantity
    - A `required` Euler angle, :math:`\alpha`. For example, ``2.1 rad``.

  * - ``beta``
    - ScalarQuantity
    - A `required` Euler angle, :math:`\beta`. For example, ``23.5 deg``.

  * - ``gamma``
    - ScalarQuantity
    - An `optional` Euler angle, :math:`\gamma`. For example, ``0.5 rad``.
      The default value is ``0 rad``.

**Example**

Using python `dict <https://docs.python.org/3/library/stdtypes.html?highlight=dict#dict>`_
object, Orientation may be represented as,

.. code-block:: py

    orientation_example = dict(alpha="0.5 rad", beta="0.23 rad", gamma="2.54 rad")





.. _symmetric_tensor:

SymmetricTensor
---------------

We define `SymmetricTensor` as an object with attributes representing the
traceless second-rank symmetric irreducible interaction tensor.

.. cssclass:: table-bordered table-hover

.. list-table:: Attributes of a SymmetricTensor object.
  :widths: 15 20 65
  :header-rows: 1

  * - Attribute name
    - Type
    - Description

  * - ``zeta``

      or

      ``Cq``

    - ScalarQuantity
    - A `required` quantity.

      **Nuclear shielding:** The strength of the anisotropy, ``zeta``, calculated
      using Haeberlen convention. The value is a physical quantity given in
      dimensionless frequency ratio, for example, ``10 ppm`` or ``10 µHz/Hz``.

      **Electric quadrupole:** The quadrupole coupling constant, ``Cq``. The
      value is a physical quantity given in units of frequency, for example,
      ``3.1 MHz``.

  * - ``eta``
    - Float
    - A `required` asymmetry parameter calculated using Haeberlen convention, for
      example, ``0.75``.

  * - ``orientation``
    - :ref:`orientation`
    - An `optional` object. The default value is a ``NULL`` object.

**Example**

Using python `dict <https://docs.python.org/3/library/stdtypes.html?highlight=dict#dict>`_
object, SymmetricTensor may be represented as,

.. code-block:: py

    # For nuclear shielding tensor
    symmetric_shielding_tensor_example = dict(
        zeta="10.3 ppm", eta=0.5, orientation=orientation_example
    )

    # For electric quadrupole tensor
    quadrupole_tensor_example = dict(
        Cq="1.45 MHz", eta=0.13, orientation=orientation_example
    )

where ``orientation_example`` is the dict object with Euler angles from the
previous example.



.. _site:

Site
----

We define a `Site` as an object with attributes representing an isolated
nuclear site.

.. cssclass:: table-bordered table-hover
.. list-table::  Attributes of a Site object.
  :widths: 25 20 55
  :header-rows: 1

  * - Attribute name
    - Type
    - Description

  * - ``isotope``
    - String
    - A `required` isotope string given as the atomic number followed by
      the isotope symbol, for example, ``13C``, ``29Si``, ``27Al``, and so on.

  * - ``isotropic_chemical_shift``
    - ScalarQuantity
    - An `optional` physical quantity describing the isotropic chemical shift
      of the site. The value is given in dimensionless frequency ratio,
      for example, ``10 ppm`` or ``10 µHz/Hz``. The default value is ``0 ppm``.

  * - ``shielding_symmetric``
    - :ref:`symmetric_tensor`
    - An `optional` object. The default is a ``NULL`` object.
      See the description for the :ref:`symmetric_tensor` object.

  * - ``quadrupolar``
    - :ref:`symmetric_tensor`
    - An `optional` object. The default is a ``NULL`` object.
      See the description for the :ref:`symmetric_tensor` object.

**Example**

Using python `dict <https://docs.python.org/3/library/stdtypes.html?highlight=dict#dict>`_
object, Site may be represented as,

.. code-block:: py

  site_example1 = dict(
      isotope="27Al",
      isotropic_chemical_shift="15 ppm",
      quadrupolar=quadrupole_tensor_example
  )

  site_example2 = dict(
      isotope="13C",
      isotropic_chemical_shift="15 ppm",
      shielding_symmetric=symmetric_shielding_tensor_example
  )

where ``symmetric_shielding_tensor_example`` and ``quadrupole_tensor_example``
are the dict object with symmetric tensor attributes from the previous example.


.. _isotopomer:

Isotopomer
----------

We define an `Isotopomer` as an object with attributes representing an isolated
spin-system.

.. An `Isotopomer` object is a python
.. `dict <https://docs.python.org/3/library/stdtypes.html?highlight=dict#dict>`__
.. object which represents an isotopomer.
.. In `mrsimulator`, each `isotopomer` is treated as a :math:`n`-coupled spin
.. system where :math:`n` is the number of sites in the isotopomer.
.. It is recommended that if the sites are uncoupled, it be specified as
.. individual isotopomers with a single site object, rather than a single
.. isotopomer with multiple sites.

.. The key-value pairs of the `Isotopomer` object follow,

.. cssclass:: table-bordered table-hover

.. list-table:: Attributes of an Isotopomer object.
  :widths: 15 15 70
  :header-rows: 1

  * - Attributes
    - Type
    - Description

  * - ``name``
    - String
    - An `optional` attribute with a name for the isotopomer.
      The default value is an empty string.

  * - ``description``
    - String
    - An `optional` attribute with a description of the isotopomer.
      The default value is an empty string.

  * - ``sites``
    - List
    - A `required` list of :ref:`site` objects.

  * - ``couplings``
    - List
    - An `optional` list of coupling objects. The default value is an empty list.
      Not yet implemented.

  * - ``abundance``
    - String
    - An `optional` quantity representing the abundance of the isotopomer. The
      abundance is given as percent, for example, ``25.4 %``. This value is useful
      when multiple isotopomers are present. The default value is ``100 %``.

**Example**

Using python `dict <https://docs.python.org/3/library/stdtypes.html?highlight=dict#dict>`_
object, Isotopomer may be represented as,

.. code-block:: py

  isotopomer_example1 = dict(sites=[site_example1], abundance="15.3 %")
  isotopomer_example2 = dict(sites=[site_example2], abundance="65.19 %")

where ``site_example1`` and ``site_example2`` are the dict objects from the
previous example.
