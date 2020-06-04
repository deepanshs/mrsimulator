
Other Objects
=============

.. _symmetric_tensor_api:

SymmetricTensor
---------------

.. currentmodule:: mrsimulator.tensors

.. autoclass:: SymmetricTensor
    :show-inheritance:

    .. rubric:: Method Documentation
    .. automethod:: to_freq_dict
    .. automethod:: to_dict_with_units


.. _antisymmetric_tensor_api:

AntisymmetricTensor
-------------------

.. currentmodule:: mrsimulator.tensors

.. autoclass:: AntisymmetricTensor
    :show-inheritance:

    .. rubric:: Method Documentation
    .. automethod:: to_freq_dict
    .. automethod:: to_dict_with_units


.. _isotope_api:

Isotope
-------

.. currentmodule:: mrsimulator.isotope

.. autoclass:: Isotope
    :show-inheritance:

    .. rubric:: Attribute Description

    .. autosummary::
        spin
        natural_abundance
        gyromagnetic_ratio
        atomic_number
        quadrupole_moment

    .. rubric:: Method Documentation

    .. automethod:: to_dict_with_units

    Example:

    .. doctest::

        >>> # 13C isotope information
        >>> carbon = Isotope(symbol='13C')
        >>> carbon.spin
        0.5
        >>> carbon.natural_abundance # in %
        1.11
        >>> carbon.gyromagnetic_ratio # in MHz/T
        10.7084
        >>> carbon.atomic_number
        6
        >>> carbon.quadrupole_moment # in eB
        0.0
