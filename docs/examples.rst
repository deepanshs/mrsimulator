

========
Examples
========

-------------------------------
Chemical Shift Anisotropy (CSA)
-------------------------------

The following examples demonstrates the simulation of chemical shift Anisotropy
spectrum. Import the `Simulator` class and the `one_D_spectrum` method from
the `mrsimulator` package as follows

.. doctest::

    >>> from mrsimulator import Simulator
    >>> from mrsimulator.methods import one_D_spectrum

To setup a simulation create an instance of the ``Simulator`` class with

.. doctest::

    >>> sim1 = Simulator()

The `sim1` instance has two attributes, isotopomers and spectrum. By default
value of isotopomers is an empty list,

    >>> print(sim1.isotopomers)
    []
and the value of the spectrum attribute is an empty dictionary.

    >>> print(sim1.spectrum)

requires two parameters to generate the spectrum,
**a list of isotopomers**, and the **spectrum metadata dictionary**.

The following in a list of isotopomers with two isotopomers, each with a
single site representing a :math:`^{13}\mathrm{C}` nucleus.

.. doctest::

    >>> isotopomers = [
    ...     {
    ...         "sites": [
    ...             {
    ...                 "isotope_symbol": "13C",
    ...                 "isotropic_chemical_shift": "1 Hz",
    ...                 "shielding_symmetric": {
    ...                     "anisotropy": "-3.89 kHz",
    ...                     "asymmetry": 0.25
    ...                 }
    ...             }
    ...         ]
    ...     },
    ...     {
    ...         "sites": [
    ...             {
    ...                 "isotope_symbol": "13C",
    ...                 "isotropic_chemical_shift": "1 kHz",
    ...                 "shielding_symmetric": {
    ...                     "anisotropy": "8.2 kHz",
    ...                     "asymmetry": 0.0
    ...                 }
    ...             }
    ...         ]
    ...     }
    ... ]

To add the list of isotopomers, use the :meth:`~mrsimulator.Simulator.set_isotopomers`
method as follows,

.. doctest::

    >>> sim1.set_isotopomers(isotopomers)



        {
            "sites": [
                {
                    "isotope_symbol": "13C",
                    "isotropic_chemical_shift": "1 kHz",
                    "shielding_symmetric": {
                        "anisotropy": "8.2 kHz",
                        "asymmetry": 0.0
                    }
                }
            ]
        },
        {
            "sites": [
                {
                    "isotope_symbol": "1H",
                    "isotropic_chemical_shift": "3 kHz",
                    "shielding_symmetric": {
                        "anisotropy": "23.2 kHz",
                        "asymmetry": 0.0
                    }
                }
            ]
        },
        {
            "sites": [
                {
                    "isotope_symbol": "29Si",
                    "isotropic_chemical_shift": "1.64 kHz",
                    "shielding_symmetric": {
                        "anisotropy": "7.36 kHz",
                        "asymmetry": 0.0
                    }
                }
            ]
        },
        {
            "sites": [
                {
                    "isotope_symbol": "29Si",
                    "isotropic_chemical_shift": "43 kHz",
                    "shielding_symmetric": {
                        "anisotropy": "8.36 kHz",
                        "asymmetry": 0.5
                    }
                }
            ]
        },
        {
            "sites": [
                {
                    "isotope_symbol": "29Si",
                    "isotropic_chemical_shift": "10 kHz",
                    "shielding_symmetric": {
                        "anisotropy": "6.36 kHz",
                        "asymmetry": 0.0
                    }
                }
            ]
        },
        {
            "sites": [
                {
                    "isotope_symbol": "1H",
                    "isotropic_chemical_shift": "5.6 kHz",
                    "shielding_symmetric": {
                        "anisotropy": "13.2 kHz",
                        "asymmetry": 0.0
                    }
                }
            ]
        }
    ]
}
