
.. _examples_from_exp:

Examples
========

In this section, we use the tools we learned so far to create isotopomers
with practical/experimental applications.

.. doctest::

    >>> from mrsimulator import Simulator, Isotopomer, Site, Dimension
    >>> from mrsimulator import SymmetricTensor as st
    >>> from mrsimulator.methods import one_d_spectrum

.. _example_coesite:
.. include:: coesite.rst

.. _example_wollastonite:
.. include:: wollastonite.rst

.. _example_potassium_sulfate:
.. include:: potassium_sulfate.rst

.. only:: html

    .. bibliography:: ref.bib
        :style: plain
