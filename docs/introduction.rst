
.. _isotopomers_docs:

********************************************
Introduction to Isotopomers and Spin systems
********************************************

Malcolm H. Levitt defines isotopomers, in his book “Spin Dynamics,” as
“Molecules differing only in the mass numbers of the nuclei are called
isotopomer.”
We can, however, generalize the concept of isotopomers by replacing the term
`molecules` with `structural units` as “Structural units differing only in the
mass numbers of the nuclei.”

It is best to illustrate isotopomers using examples. Consider a structural
unit, H-C. The most abundant isotopes of H and C are :math:`^1\text{H}`
(99.985%), :math:`^2\text{H}` (0.015%), and :math:`^{12}\text{C}` (98.93%),
:math:`^{13}\text{C}` (1.11%), respectively, where we consider only the top two
most abundant isotopes. From this, we can create four H-C isotopomers, as
listed in :numref:`isotopomers_list`. Here, each isotopomer consists of two
isotopes. From an NMR viewpoint, the isotopomers at index 1 and 2, are
considered as a single-site spin system, because the corresponding isotope of
carbon, :math:`^{12}\text{C}`, is NMR inactive. The
isotopomers at index 3 and 4 are two-site spin system with a single coupling.
Here, the isotopomer abundance is given as the product of the natural abundance of
the individual isotopes, compositing the isotopomer.
The observed NMR signal is the sum of the signals arising from individual
isotopomers, weighted by their respective abundance.

.. cssclass:: table-bordered table-striped centered
.. _isotopomers_list:
.. list-table:: Four isotopomers resulting from H-C structural unit.
   :widths: 15 15 15 15 40
   :header-rows: 1

   * - Index
     - Isotopomers
     - Sites
     - Coupled pairs
     - Abundance

   * - 1
     - 1H - 12C
     - 1
     - 0
     - (99.985% x 98.93%) ~ 98.915%

   * - 2
     - 2H - 12C
     - 1
     - 0
     - (0.015% x 98.93%) ~ 0.148%

   * - 3
     - 1H - 13C
     - 2
     - 1
     - (99.985% x 1.11%) ~ 1.11%

   * - 4
     - 2H - 13C
     - 2
     - 1
     - (0.015% x 1.11%) ~ 0.00016%

In the ``mrsimulator`` library, we consider each isotopomer as a simplified isolated
spin system, where only the NMR active sites may reside within the spin system.
All NMR inactive sites are ignored. This simplified isolated spin system is given
the class name of **SpinSystem**. In the following sub-section, we illustrate with
examples of how we represent the isotopomers using the spin systems. For a detailed
description of the class attributes, refer to :numref:`table_spin_system` to
:numref:`table_symmetric_tensor`.


Overview of the SpinSystem Model
--------------------------------

Uncoupled spin systems
''''''''''''''''''''''

.. _listing_1H-12C:
.. code-block:: json
  :linenos:
  :caption: An example 1H-12C isotopomer in JSON representation.

  {
      "name": "1H-12C",
      "description": "An optional description of the spin system/isotopomer",
      "sites": [
          {
              "isotope": "1H",
              "isotropic_chemical_shift": "-1.2 ppm",
              "shielding_symmetric": {
                  "zeta": "4.12 ppm",
                  "eta": 0.12
              }
          }
      ],
      "couplings": [],
      "abundance": "98.915%"
  }

:numref:`listing_1H-12C` is an example of the spin system corresponding to the
`1H-12C` isotopomer, serialized using the JavaScript Object Notation (JSON). At the
root level of the **SpinSystem** object, we find five keywords, **name**,
**description**, **sites**, **couplings**, and **abundance**. The value of the `name`
key is the name of the spin system/isotopomer, here given a value of `1H-12C`. The
value of the description key is an optional string describing the spin system. The
value of the `sites` key is a list of **Site** objects. Here, this list comprises of a
single **Site** object (lines 5-12).
The value of the `couplings` key is a list of **Coupling** objects. In this example,
there are no  couplings, and hence the value of this attribute is an empty list.
The value of the `abundance` key is the abundance of the spin system, here given a
value of `98.915%` based on the data from :numref:`isotopomers_list`.
See :numref:`table_spin_system` for further description of the **SpinSystem** class and
its attributes.

The **Site** object (lines 5-12) is described with three keywords, **isotope**,
**isotropic_chemical_shift**, and **shielding_symmetric**. The value of the `isotope`
key is the spin isotope, here given a value of `1H`.
The value of the `isotropic_chemical_shift`, `-1.2 ppm`, is the
:math:`^1\text{H}` isotropic chemical shift. Because :math:`^1\text{H}` is
:math:`I = 1/2`, we have additionally defined an optional `shielding_symmetric`,
which represents the components of the second-rank traceless symmetric nuclear
shielding tensor. We parameterize this tensor using the Haeberlen convention with
parameters `zeta` and `eta`, defined as the strength of the anisotropy and asymmetry,
respectively. See :numref:`table_site` for further information on the **Site** object
and its attributes.


.. _listing_2H-12C:
.. code-block:: json
  :linenos:
  :emphasize-lines: 12-18
  :caption: An example 2H-12C isotopomer in JSON representation.

  {
      "name": "2H-12C",
      "description": "An optional description on the isotopomer",
      "sites": [
          {
              "isotope": "2H",
              "isotropic_chemical_shift": "4.1 ppm",
              "shielding_symmetric": {
                  "zeta": "12.12 ppm",
                  "eta": 0.82
              },
              "quadrupolar": {
                  "Cq": "1.47 MHz",
                  "eta": 0.27,
                  "alpha": "0.212 rad",
                  "beta": "1.231 rad",
                  "gamma": "3.1415 rad"
              }
          }
      ],
      "coupling": [],
      "abundance": "0.148%"
  }

:numref:`listing_2H-12C` is an example of a spin system representing the `2H-12C`
isotopomer. This example is similar to the example from :numref:`listing_1H-12C`,
except we see a new keyword, **quadrupolar**, in the **Site** object (lines 12-18).
In this example, the site `isotope` is `2H`, which is a quadrupolar nucleus,
:math:`I>1/2`. For quadrupolar nuclei, besides the nuclear shielding tensor, there
also exists an electric field gradient (EFG) tensor. An EFG tensor is a second-rank
traceless symmetric tensor, which we describe by the parameters `Cq` and
`eta` as the quadrupolar coupling constant and asymmetry parameter, respectively.
Additionally, we see the Euler angle orientations, `alpha`, `beta`, and `gamma`, which
are the relative orientation of the EFG tensor from the nuclear shielding tensor.


Coupled spin systems
''''''''''''''''''''

.. note::
    The current version of the ``mrsimulator`` package does not include coupled
    spin systems. The SpinSystem model for the couplings will be made available when
    we include the coupled spin systems to the package.


Table of Class Attributes
-------------------------

.. cssclass:: table-bordered table-striped centered
.. _table_spin_system:
.. list-table:: The attributes of a SpinSystem object.
  :widths: 15 15 70
  :header-rows: 1

  * - Attributes
    - Type
    - Description

  * - ``name``
    - String
    - An `optional` attribute with a name for the isotopomer/spin system. Naming is a
      good practice as it improves the readability, especially when multiple
      spin systems are present. The default value is an empty string.

  * - ``description``
    - String
    - An `optional` attribute describing the spin system. The default value is an empty
      string.

  * - ``sites``
    - List
    - An `options` list of :ref:`site` objects. The default value is an empty list.

  * - ``couplings``
    - List
    - An `optional` list of coupling objects. The default value is an empty list.
      Not yet implemented.

  * - ``abundance``
    - String
    - An `optional` quantity representing the abundance of the isotopomer/spin system.
      The abundance is given as percentage, for example, ``25.4 %``. This value is
      useful when multiple spin systems are present. The default value is ``100 %``.


.. cssclass:: table-bordered table-striped centered
.. _table_site:
.. list-table::  The attributes of a Site object.
  :widths: 30 15 50
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
    - An `optional` object describing the second-rank traceless symmetric
      nuclear shielding tensor following the Haeberlen convention. The default is a
      ``NULL`` object. See the description for the :ref:`symmetric_tensor` object.

  * - ``quadrupolar``
    - :ref:`symmetric_tensor`
    - An `optional` object describing the second-rank traceless electric
      quadrupole tensor. The default is a ``NULL`` object.
      See the description for the :ref:`symmetric_tensor` object.



.. cssclass:: table-bordered table-striped centered
.. _table_symmetric_tensor:
.. list-table:: The attributes of a SymmetricTensor object.
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
      using the Haeberlen convention. The value is a physical quantity given in
      dimensionless frequency ratio, for example, ``10 ppm`` or ``10 µHz/Hz``.

      **Electric quadrupole:** The quadrupole coupling constant, ``Cq``. The
      value is a physical quantity given in units of frequency, for example,
      ``3.1 MHz``.

  * - ``eta``
    - Float
    - A `required` asymmetry parameter calculated using the Haeberlen convention, for
      example, ``0.75``.

  * - ``alpha``
    - ScalarQuantity
    - An `optional` Euler angle, :math:`\alpha`. For example, ``2.1 rad``.
      The default value is ``0 rad``.

  * - ``beta``
    - ScalarQuantity
    - An `optional` Euler angle, :math:`\beta`. For example, ``90°``.
      The default value is ``0 rad``.

  * - ``gamma``
    - ScalarQuantity
    - An `optional` Euler angle, :math:`\gamma`. For example, ``0.5 rad``.
      The default value is ``0 rad``.
