


.. _isotopomers:

------------------
Isotopomers object
------------------

An `Isotopomers` object is a python list object with a collection of
**Isotopomer** objects.



**Isotopomer object**
---
  An `Isotopomer` object is a python dict object with the following
  key-value pairs,

  .. list-table::
    :widths: 25 75
    :header-rows: 1

    * - key
      - value description
    * - ``sites``
      - The value is a list of **Site** objects.
    * - ``abundance``
      - The value is a float with the abundance of the isotopomer.
        This is useful where more than one isotopomers are simulated.

  ..  * - ``coulpings``
  ..    - Not yet implemented.


  In `mrsimulator` each isotopomer is treated as an :math:`n`-coupled spin
  problem where :math:`n` is the number of sites in the isotopomer
  It is recomended that if the sites are uncoupled, they be specified as
  individual isotopomers rather a single isotoper with multiple sites.


**Site object**
---
  A `Site` object is again a python dict object. The key-value pairs of this
  object describes a nuclear site.

  .. list-table::
    :widths: 25 75
    :header-rows: 1

    * - key
      - value description
    * - ``isotope_symbol``
      - The value is a string representing the NMR active isotope
        symbol. In the above example, this value is '13C'.
    * - ``isotropic_chemical_shift``
      - The value is a string containing a physical quantity that represents
        the isotropic chemical shift of the NMR active isotope symbol, for example,
        "1 Hz".
    * - ``shielding_symmetric``
      - The value is a `SymmetricTensor` object describing the parameters of a
        traceless second rank symmetric tensor.



**SymmetricTensor object**




**Example of Isotopomers object**

As an example, consider the following `Isotopomers` object,

    >>> [
    ...     {
    ...         "sites": [
    ...             {
    ...                 "isotope_symbol": "13C",
    ...                 "isotropic_chemical_shift": "1 Hz",
    ...                 "shielding_symmetric": {
    ...                     "anisotropy": "-3.89 kHz",
    ...                     "asymmetry": 0.
    ...                 }
    ...             }
    ...         ]
    ...     },
    ...     {
    ...         "sites": [
    ...             {
    ...                 "isotope_symbol": "29Si",
    ...                 "isotropic_chemical_shift": "0 Hz",
    ...                 "shielding_symmetric": {
    ...                     "anisotropy": "8.89 kHz",
    ...                     "asymmetry": 0.5
    ...                 }
    ...             }
    ...         ]
    ...     }
    ... ]

The above example contains two `Isotopomer` objects, each with a single `Site`
object. The first isotopomer describes a :math:`^{13}\mathrm{C}` nuclear site
with 1 Hz isotropic chemical shift, -3.89 kHz of anisotropy with an asymmetry
parameter of 0.0. The second isotopomer describes a :math:`^{29}\mathrm{Si}`
nuclear site with no isotropic shift, 8.89 kHz of anisotropy and a value of 0.5
for the asymmetry parameter.


.. Note::
    All physical quantities are specified as strings containing a numerical
    value and a unit.


Almost all interactions in nuclear magnetic resonance can be described by a
second rank tensor with nine components. This second rank tensor can be divided
into: isotropic component, the
