

.. _isotopomer_api:

Isotopomer
==========

.. currentmodule:: mrsimulator

.. autoclass:: Isotopomer
    :show-inheritance:

    .. rubric:: Method Documentation

    .. automethod:: parse_dict_with_units


.. _site_api:

Site
====

.. currentmodule:: mrsimulator

.. autoclass:: Site
    :show-inheritance:

    .. rubric:: Attribute Description
    .. autosummary::
        spin
        natural_abundance
        gyromagnetic_ratio
        atomic_number
        quadrupole_moment

    Example:

    .. doctest::

        >>> # 13C isotope information
        >>> site1.spin
        0.5
        >>> site1.natural_abundance
        1.11
        >>> site1.gyromagnetic_ratio
        10.7084
        >>> site1.atomic_number
        6
        >>> site1.quadrupole_moment
        0.0

.. _symmetric_tensor_api:

SymmetricTensor
===============

.. currentmodule:: mrsimulator

.. autoclass:: SymmetricTensor
    :show-inheritance:

    .. rubric:: Method Documentation
    .. automethod:: to_freq_dict
    .. automethod:: to_dict_with_units


.. _antisymmetric_tensor_api:

AntisymmetricTensor
===================

.. currentmodule:: mrsimulator

.. autoclass:: AntisymmetricTensor
    :show-inheritance:

    .. rubric:: Method Documentation
    .. automethod:: to_freq_dict
    .. automethod:: to_dict_with_units


.. _dimension_api:

Dimension
=========

.. currentmodule:: mrsimulator

.. autoclass:: Dimension
    :show-inheritance:

    .. autoattribute:: coordinates_Hz
    .. autoattribute:: coordinates_ppm


    .. rubric:: Method Documentation

    .. automethod:: parse_dict_with_units
