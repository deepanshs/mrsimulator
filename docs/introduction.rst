
.. _isotopomers_docs:

****************************
Introduction to Spin Systems
****************************


At the heart of any ``mrsimulator`` calculation is the definition of a **SpinSystem** object describing the sites and couplings within a spin system.  We begin by examining the definition of a **Site** object. 

Site
''''

Consider the example below of the JSON serialization of a **Site** object for a deuterium nucleus.

.. _listing_1H:
.. code-block:: json
  :linenos:
  :caption: An example 2H site in JSON representation.

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

The value of the `isotope` key holds the spin isotope, here given a value of `2H`. 
The value of the `isotropic_chemical_shift` is the optional :math:`^2\text{H}` isotropic chemical shift, here given as `4.1 ppm`. We have additionally defined an optional `shielding_symmetric`, whose value holds a dictionary with the components of the second-rank traceless symmetric nuclear
shielding tensor. We parameterize this tensor using the Haeberlen convention with
parameters `zeta` and `eta`, defined as the strength of the anisotropy and asymmetry,
respectively.  Since deuterium is a quadrupolar nucleus, :math:`I>1/2`, there also can
be a quadrupolar coupling interaction between the nuclear quadrupole moment and the
surrounding electric field gradient (EFG) tensor, defined in a dictionary held in the optional 
key `quadrupolar`.  An EFG tensor is a second-rank traceless symmetric tensor, and we describe
the quadrupolar coupling with the parameters `Cq` and `eta`, i.e., the quadrupolar coupling constant 
and asymmetry parameter, respectively.  Additionally, we see the Euler angle orientations, 
`alpha`, `beta`, and `gamma`, which are the relative orientation of the EFG tensor from the 
nuclear shielding tensor. 


See :numref:`table_site` and :numref:`table_symmetric_tensor` for further information on the **Site** and **SymmetricTensor** objects and their attributes, respectively.


Table of Site Class Attributes
------------------------------

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

SpinSystem
''''''''''

As mentioned earlier, the **SpinSystem** object, used in the ``mrsimulator`` package, describes the sites and couplings within a spin system. 


Uncoupled spin systems
----------------------

Using the previous **Site** object example, we construct a simple single site **SpinSystem** object shown below.

.. _listing_2H:
.. code-block:: json
  :linenos:
  :caption: An example 2H spin system in JSON representation.

  {
      "name": "2H spin system",
      "description": "An optional description on the spin system",
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
      "couplings": [],
      "abundance": "0.148%"
  }
  
At the root level of the **SpinSystem** object, we find four keywords, **name**,
**description**, **sites**, and **abundance**. The value of the `name`
key is the name of the spin system, here given a value of `2H spin system`. The
value of the description key is an optional string describing the spin system. The
value of the `sites` key is a list of **Site** objects. Here, this list comprises of a
single **Site** object (lines 5-19). The value of the `abundance` key is the 
abundance of the spin system, here given a value of `0.148%`.
The value of the `couplings` key is a list 
of **Coupling** objects. In this example, there are no  couplings, and hence the value of 
this attribute is an empty list. See :numref:`table_spin_system` for further 
description of the **SpinSystem** class and its attributes.




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
    - An `optional` attribute with a name for the spin system. Naming is a
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
    - An `optional` quantity representing the abundance of the spin system.
      The abundance is given as percentage, for example, ``25.4 %``. This value is
      useful when multiple spin systems are present. The default value is ``100 %``.


Coupled spin systems
----------------------

.. note::
    The current version of the ``mrsimulator`` package does not include coupled
    spin systems. The SpinSystem model for the couplings will be made available when
    we include the coupled spin systems to the package.  The ``mrsimulator`` package
    will eventually handle coupled spin systems, but only in the weak coupling limit.


