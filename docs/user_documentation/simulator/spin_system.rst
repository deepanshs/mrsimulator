.. _spin_system_documentation:

**********************
Spin Systems and Sites
**********************


At the heart of any ``mrsimulator`` calculation is the definition of a **SpinSystem**
object describing the sites and couplings within a spin system. Each **Simulator** holds
a list of **SpinSystem** objects which will be used to calculate the spectrum. We begin
by examining the definition of a **Site** object.

Site
''''

A site object holds single site NMR interaction parameters, which include the nuclear
shielding and quadrupolar interaction parameters.
Consider the example below of a **Site** object for a deuterium nucleus created in Python.

.. _listing_site:
.. code-block:: python

    # Import objects for the Site
    from mrsimulator import Site
    from mrsimulator.spin_system.tensors import SymmetricTensor

    # Create the site object
    H2_site = Site(
        isotope="2H",
        isotropic_chemical_shift=4.1,  # in ppm
        shielding_symmetric=SymmetricTensor(
            zeta=12.12,  # in ppm
            eta=0.82,
            alpha=5.45,  # in radians
            beta=4.82,   # in radians
            gamma=0.5,   # in radians
        ),
        quadrupolar=SymmetricTensor(
            Cq=1.47e6,     # in Hz
            eta=0.27,
            alpha=0.212,   # in radians
            beta=1.231,    # in radians
            gamma=3.1415,  # in radians
        )
    )

The *isotope* key holds the spin isotope, here given a value of *2H*.
The *isotropic_chemical_shift* is the isotropic chemical shift of the site isotope, :math:`^2\text{H}`,
here given as *4.1 ppm*. We have additionally defined an optional *shielding_symmetric* key, whose
value is a second-rank traceless symmetric nuclear shielding
tensor.

.. note::
  We parameterize a SymmetricTensor using the Haeberlen convention with parameters *zeta* and *eta*,
  defined as the shielding anisotropy and asymmetry, respectively. The Euler angle orientations, *alpha*,
  *beta*, and *gamma* are the relative orientation of the nuclear shielding tensor from a common reference
  frame.

Since deuterium is a quadrupolar nucleus, :math:`I>1/2`, there also can be a quadrupolar coupling
interaction between the nuclear quadrupole moment and the surrounding electric field gradient (EFG) tensor,
defined in the optional *quadrupolar* key. An EFG tensor is a second-rank traceless
symmetric tensor, and we describe its components with *Cq* and *eta*, i.e., the quadrupolar coupling constant
and asymmetry parameter, respectively.  Additionally, we use the Euler angle orientations, *alpha*, *beta*,
and *gamma*, which are the relative orientation of the EFG tensor from a common reference frame.

See :numref:`table_site` and :numref:`table_symmetric_tensor` for further information on
the **Site** and **SymmetricTensor** objects and their attributes, respectively. Also, all ``mrsimulator``
objects have the attribute *property_units* which provides the units for all class properties

.. code_block:: python

    print(Site().property_units)
    # {'isotropic_chemical_shift': 'ppm'}

    print(SymmetricTensor().property_units)
    # {'zeta': 'ppm', 'Cq': 'Hz', 'D': 'Hz', 'alpha': 'rad', 'beta': 'rad', 'gamma': 'rad'}


SpinSystem
''''''''''

The **SpinSystem** object is a collection of sites and couplings within a spin system.

Single Site Spin System
-----------------------

Using the previous 2H **Site** object example, we construct a simple single-site
**SpinSystem** object, as shown below.

.. _listing_single_site_spin_system:
.. code-block:: python

    # Import the Spin System object
    from mrsimulator import SpinSystem

    single_site_sys = SpinSystem(
        name="2H uncoupled spin system",
        description="An optional description of the spin system",
        sites=[H2_site],
        abundance=80,   # percentage
    )

At the root level of the **SpinSystem** object, we find four keywords, **name**,
**description**, **sites**, and **abundance**. The value of the *name* key is the
optional name of the spin system. The
value of the description key is an optional string describing the spin system. The
value of the *sites* key is a list of **Site** objects. Here, this list is simply
the single single **Site** object, `2H_site`, from :ref:`the Site object <listing_site>`.
The value of the *abundance* key is the abundance of the spin system, here given
a value of *80%*.

See :numref:`table_spin_system` for further description of the **SpinSystem** class and
its attributes.

Multi Site Spin System
----------------------

To create a **SpinSystem** with more than one **Site**, the *sites* attribute simply needs
to be set to a list of **Site** objects. Here we create a :math:`^{13}\text{C}` site and add
it to the *sites* attribute of our already defined **SpinSystem**.

.. _listing_multi_site_spin_system:
.. code-block:: python
    :format: doctest
    :context: close-figs
    :include-source:

    # Create the new Site object
    C13_site = Site(
        isotope="13C",
        isotropic_chemical_shift=-53.2,  # in ppm
        shielding_symmetric=SymmetricTensor(
                zeta=90.5,  # in ppm
                eta=0.64,
            )
        )
    )

    # Create a new SpinSystem object with both Sites
    multi_site_sys = SpinSystem(
        name="Multi site spin system",
        description="A spin system with multiple sites",
        sites=[H2_site, C13_site],
        abundance=0.148,  # percentage
    )

Again we see the optional *name* and *description* attributes. The *sites* attribute is now
a list of two **Site** objects, the previous :math:`^2\text{H}` site and the new
:math:`^{13}\text{C}` site.

Coupling
''''''''

A coupling object holds two site NMR interaction parameters, which includes the *J*-coupling
and the dipolar coupling interaction parameters.
Consider the example below of a **Coupling** object between the two sites in for the
previous **SpinSystem**.

.. _listing_coupling:
.. code-block:: python

    # Import the Coupling object
    from mrsimulator import Coupling

    H2_C13_coupling = Coupling(
        site_index=[0, 1],
        isotropic_j=15,  # in Hz
        j_symmetric=SymmetricTensor(
            zeta=12.12,  # in Hz
            eta=0.82,
            alpha=2.45,  # in radians
            beta=1.75,   # in radians
            gamma=0.15,  # in radians
        ),
        dipolar=SymmetricTensor(
            D=1.7e3,       # in Hz
            alpha=0.12,    # in radians
            beta=0.231,    # in radians
            gamma=1.1415,  # in radians
        )
    )

The *site_index* key holds a list of two integers corresponding to the index of the two coupled sites
within the spin system. We will add this coupling object to ``multi_site_sys`` so 0 and 1 correspond to
the :math:`^2\text{H}` and :math:`^{13}\text{C}`  The value of the *isotropic_j* is the isotropic
*J*-coupling, here given as *15 Hz*. We have additionally defined an optional *j_symmetric* key,
whose value holds a SymmetricTensor object representing the traceless 2nd-rank symmetric *J*-coupling
tensor.

Additionally, the dipolar coupling interaction between the coupled nuclei is defined with an optional
*dipolar* key. A dipolar tensor is a second-rank traceless symmetric tensor, and we describe the dipolar
coupling constant with the parameter *D*.  The Euler angle orientations, *alpha*, *beta*, and *gamma*
are the relative orientation of the dipolar tensor from a common reference frame.

See :numref:`table_coupling` and :numref:`table_symmetric_tensor` for further information on
the **Site** and **SymmetricTensor** objects and their attributes, respectively.

Coupled Spin System
-------------------

By adding the coupling object to ``multi_site_sys``, we create a coupled spin system.

.. _listing_coupled_spin_system:
.. code-block:: python

    # We have no couplings present in the SpinSystem
    print(multi_site_sys.couplings)
    # []

    # Set the coupling attribute to the previously defined Coupling object
    multi_site_sys.couplings = [H2_13C_coupling]

In comparison to the previous example, there are five keywords at the root level of the
**SpinSystem** object, **name**, **description**, **sites**, **couplings**, and **abundance**.
In this example, the value of the *sites* key holds a list of two **Site** objects. At index 0
(lines 9-26) is the *2H* site from :ref:`the Site object <listing_site>`, and at index 1
(lines 27-34) is a *13C* site. The value of the *couplings* key is a list of **Coupling** objects.
In this example, we define a single coupling object (lines 37-60) from the
:ref:`coupling object <listing_coupling>`. The value of the *site_index* key within the
coupling object, line 34, corresponds to the site index from the *sites* list.


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
    - An *optional* attribute with a name for the spin system. Naming is a
      good practice as it improves the readability, especially when multiple
      spin systems are present. The default value is an empty string.

  * - ``label``
    - String
    - An *optional* attribute giving a label to the spin system. Like ``name``, it has no
      effect on a simulation and is purely for readability.

  * - ``description``
    - String
    - An *optional* attribute describing the spin system. The default value is an empty
      string.

  * - ``sites``
    - List
    - An *optional* list of :ref:`site_api` objects. The default value is an empty list.

  * - ``couplings``
    - List
    - An *optional* list of coupling objects. The default value is an empty list.

  * - ``abundance``
    - String
    - An *optional* quantity representing the abundance of the spin system.
      The abundance is given as percentage, for example, ``25.4`` for 25.4%. This value is
      useful when multiple spin systems are present. The default value is ``100``.


.. cssclass:: table-bordered table-striped centered
.. _table_site:
.. list-table::  The attributes of a Site object.
  :widths: 30 15 50
  :header-rows: 1

  * - Attribute name
    - Type
    - Description

  * - ``name``, ``label``, and ``description``
    - String
    - All three are *optional* attributes giving context to a **Site** object. The default
      value for all three is an empty string.

  * - ``isotope``
    - String
    - A *required* isotope string given as the atomic number followed by
      the isotope symbol, for example, ``13C``, ``29Si``, ``27Al``, and so on.

  * - ``isotropic_chemical_shift``
    - ScalarQuantity
    - An *optional* physical quantity describing the isotropic chemical shift
      of the site. The value is given in ppm, for example, ``10`` for 10 ppm.
      The default value is ``0``.

  * - ``shielding_symmetric``
    - :ref:`sy_api`
    - An *optional* object describing the second-rank traceless symmetric
      nuclear shielding tensor following the Haeberlen convention. The default
      is ``None``. See the description for the :ref:`sy_api` object.

  * - ``quadrupolar``
    - :ref:`sy_api`
    - An *optional* object describing the second-rank traceless electric
      quadrupole tensor. The default is ``None``.
      See the description for the :ref:`sy_api` object.


.. cssclass:: table-bordered table-striped centered
.. _table_coupling:
.. list-table::  The attributes of a Coupling object.
  :widths: 30 15 50
  :header-rows: 1

  * - Attribute name
    - Type
    - Description

  * - ``site_index``
    - List of two integers
    - A *required* list with integers corresponding to the site index of the coupled
      sites, for example, [0, 1], [2, 1]. The order of the integers is irrelevant.

  * - ``isotropic_j``
    - ScalarQuantity
    - An *optional* physical quantity describing the isotropic *J*-coupling in Hz.
      The default value is ``0``.

  * - ``j_symmetric``
    - :ref:`sy_api`
    - An *optional* object describing the second-rank traceless symmetric *J*-coupling
      tensor following the Haeberlen convention. The default is ``None``. See
      the description for the :ref:`sy_api` object.

  * - ``dipolar``
    - :ref:`sy_api`
    - An *optional* object describing the second-rank traceless dipolar tensor. The
      default is ``None``. See the description for the :ref:`sy_api`
      object.


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

      or

      ``D``

    - ScalarQuantity
    - A *required* quantity.

      **Nuclear shielding:** The shielding anisotropy, ``zeta``, calculated
      using the Haeberlen convention. The value is a physical quantity given in
      ppm, for example, ``10``

      **Electric quadrupole:** The quadrupole coupling constant, ``Cq``. The
      value is a physical quantity given in units of Hz, for example,
      ``3.1e6`` for 3.1 MHz.

      **J-coupling:** The *J*-coupling anisotropy, ``zeta``, calculated
      using the Haeberlen convention. The value is a physical quantity given in
      Hz, for example, ``10`` for 10 Hz.

      **Dipolar-coupling:** The dipolar-coupling constant, ``D``. The value is a
      physical quantity given in Hz, for example, ``9e6`` for 9 kHz.

  * - ``eta``
    - Float
    - A *required* asymmetry parameter calculated using the Haeberlen convention, for
      example, ``0.75``. The parameter is set to zero for the dipolar tensor.

  * - ``alpha``
    - ScalarQuantity
    - An *optional* Euler angle, :math:`\alpha`. For example, ``2.1`` for 2.1 radians.
      The default value is ``0``.

  * - ``beta``
    - ScalarQuantity
    - An *optional* Euler angle, :math:`\beta`. For example, ``1.5708`` for 90 degrees.
      The default value is ``0``.

  * - ``gamma``
    - ScalarQuantity
    - An *optional* Euler angle, :math:`\gamma`. For example, ``0.5`` for 0.5 radians.
      The default value is ``0``.
