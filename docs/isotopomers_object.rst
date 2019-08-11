

.. _isotopomers:

-----------
Isotopomers
-----------

The `Isotopomers` object is a python
`list <https://docs.python.org/3/library/stdtypes.html#list>`_ object
with a collection of :ref:`isotopomer` objects. The simulated lineshape
is a weighted average over all the isotopomers where the weights are the
abundance of each isotopomer.



*Example of Isotopomers object*

    >>> isotopomers_list = [
    ...     {
    ...         "sites": [
    ...             {
    ...                 "isotope": "13C",
    ...                 "isotropic_chemical_shift": "1 ppm",
    ...                 "shielding_symmetric": {
    ...                     "anisotropy": "-3.89 µHz/Hz",
    ...                     "asymmetry": 0.
    ...                 }
    ...             }
    ...         ],
    ...         "abundance": "100 %"
    ...     },
    ...     {
    ...         "sites": [
    ...             {
    ...                 "isotope": "29Si",
    ...                 "isotropic_chemical_shift": "-0.1 mHz/Hz",
    ...                 "shielding_symmetric": {
    ...                     "anisotropy": "8.89 ppm",
    ...                     "asymmetry": 0.5
    ...                 }
    ...             }
    ...         ],
    ...         "abundance": "100 %"
    ...     }
    ... ]

In the above example, the variable ``isotopomers_list`` is a list with two
:ref:`isotopomer` objects, each containing a single :ref:`site` object.
The first isotopomer describes a :math:`^{13}\mathrm{C}` nuclear site
with 1 ppm isotropic chemical shift, -3.89 µHz/Hz of anisotropy and an
asymmetry
parameter of 0.0. The second isotopomer describes a :math:`^{29}\mathrm{Si}`
nuclear site with -0.1 mHz/Hz of isotropic chemical shift, 8.89 ppm of
anisotropy
and a value of 0.5 for the asymmetry parameter.
