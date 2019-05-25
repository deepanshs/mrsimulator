
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
