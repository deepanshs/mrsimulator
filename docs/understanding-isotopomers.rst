
.. _isotopomers_docs:

***************************
Introduction to Isotopomers
***************************

Malcolm H. Levitt defines isotopomers, in his book “Spin Dynamics,” as
“Molecules differing only in the mass numbers of the nuclei are called
isotopomer.”
We can, however, generically define isotopomers by replacing the term
`molecules` with `structural units` as “Structural units differing only in the
mass numbers of the nuclei are called isotopomer.”

It is best to illustrate isotopomers using examples. Consider a structural
unit, H-C. The most abundant isotopes of H and C are :math:`^1\text{H}`
(99.985%), :math:`^2\text{H}` (0.015%), and :math:`^{12}\text{C}` (98.93%),
:math:`^{13}\text{C}` (1.11%), respectively, where we consider only the top two
most abundant isotopes. From this, we can create four H-C isotopomers, as
listed in :numref:`isotopomers_list`. Here, each isotopomer consists of two
isotopes. From an NMR viewpoint, the isotopomers at index 1 and 2, are
considered as single-site isotopomers, because the corresponding isotope of
carbon, :math:`^{12}\text{C}`, is NMR inactive. The
isotopomers at index 3 and 4 are two-site isotopomers with a single coupling.
The isotopomer abundance is given as the product of the natural abundance of
the individual isotopes, compositing the isotopomer.

The observed NMR signal is a sum of the signals arising from individual
isotopomers, weighted by the respective abundance.

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



Overview of the Isotopomer Model
--------------------------------

In designing the **Isotopomer** model, we follow a similar premise.
An isotopomer class consists of `name`, `description`, `sites`, `couplings`,
and `abundance` attributes.
This section is a general overview of the **Isotopomer** class design and its
attributes.


.. _listing_1H-12C:
.. code-block:: json
   :linenos:
   :caption: An example 1H-12C isotopomer in JSON representation.

    {
        "name": "1H-12C",
        "description": "An optional description on the isotopomer",
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
        "abundance": "98.915%"
    }

In :numref:`listing_1H-12C` is an example of the `1H-12C` isotopomer,
serialized using the JavaScript Object Notation (JSON).
At the root level of the **Isotopomer** object, we find four keywords,
**name**, **description**, **sites**, and **abundance**. The value of the
`name` key is `1H-12C`. The value of the description key is an optional
string describing the isotopomer. The value of the `sites` key is a list of
**Site** objects. Here, the list comprises of only one **Site**
object (lines 5-12). The value of the `abundance` key is the abundance of the
isotopomer, here, given a value of `98.915%` based on data from
:numref:`isotopomers_list`. See :numref:`table_isotopomer` for
further description of **Isotopomer** object and its attributes.

The **Site** object (lines 5-12) is described with three keywords, **isotope**,
**isotropic_chemical_shift**, and **shielding_symmetric**. Here, the value of
the `isotope` is `1H`. The value of the `isotropic_chemical_shift`, `-1.2 ppm`,
is the :math:`^1\text{H}` isotropic chemical shift. Because :math:`^1\text{H}`
is :math:`I = 1/2`, we have additionally defined an optional
`shielding_symmetric`,
which represents the second rank traceless symmetric nuclear shielding tensor,
using Haeberlen convention. In this example, `zeta` and `eta` are the
anisotropy strength and asymmetry parameter, respectively. See
:numref:`table_site` for further information on the **Site** object and its
attributes.


.. _listing_2H-12C:
.. code-block:: json
   :linenos:
   :emphasize-lines: 12-17
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
                    "beta": "1.231 rad"
                }
            }
        ],
        "abundance": "0.148%"
    }

In :numref:`listing_2H-12C` is an example of a `2H-12C` isotopomer. This
example is similar to the example in :numref:`listing_1H-12C`, except we have
defined a new keyword, **quadrupolar**, to the **Site** object (lines 12-17).
In this example, the site `isotope` is `2H`, which is a quadrupolar nucleus,
:math:`I>1/2`. For quadrupolar nuclei, besides nuclear shielding tensor, there
also exists an electric field gradient (EFG) tensor. An EFG tensor is a
second-rank traceless symmetric tensor, which is described here with
parameters, `Cq` and `eta`, the quadrupolar coupling constant and asymmetry
parameter, respectively. Additionally, we have also provided the Euler angle
orientation, `alpha`, and `beta`, which gives the relative orientation of the
EFG tensor with respect to the nuclear shielding tensor.


.. note::
    The current version of the `mrsimulator` package does not support coupled
    spin-systems. The isotopomer model for couplings will be made available when
    we include the coupled spin-systems to the package.


.. cssclass:: table-bordered table-striped centered
.. _table_isotopomer:
.. list-table:: The attributes of an Isotopomer object.
  :widths: 15 15 70
  :header-rows: 1

  * - Attributes
    - Type
    - Description

  * - ``name``
    - String
    - An `optional` attribute with a name for the isotopomer. Naming is a good
      practice as it improves the readability, especially when multiple
      isotopomers are defined. The default value is an empty string.

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
      nuclear shielding tensor using Haeberlen convention. The default is a
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
      using Haeberlen convention. The value is a physical quantity given in
      dimensionless frequency ratio, for example, ``10 ppm`` or ``10 µHz/Hz``.

      **Electric quadrupole:** The quadrupole coupling constant, ``Cq``. The
      value is a physical quantity given in units of frequency, for example,
      ``3.1 MHz``.

  * - ``eta``
    - Float
    - A `required` asymmetry parameter calculated using Haeberlen convention, for
      example, ``0.75``.

  * - ``alpha``
    - ScalarQuantity
    - An `optional` Euler angle, :math:`\alpha`. For example, ``2.1 rad``.
      The default value is ``0 rad``.

  * - ``beta``
    - ScalarQuantity
    - An `optional` Euler angle, :math:`\beta`. For example, ``90 deg``.
      The default value is ``0 rad``.

  * - ``gamma``
    - ScalarQuantity
    - An `optional` Euler angle, :math:`\gamma`. For example, ``0.5 rad``.
      The default value is ``0 rad``.
