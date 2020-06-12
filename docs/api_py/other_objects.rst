
Other Objects
=============

.. _symmetric_tensor_api:

SymmetricTensor
---------------

.. currentmodule:: mrsimulator.spin_system.tensors

.. autoclass:: SymmetricTensor
    :show-inheritance:

    .. rubric:: Method Documentation
    .. automethod:: to_freq_dict
    .. automethod:: to_dict_with_units


.. _antisymmetric_tensor_api:

AntisymmetricTensor
-------------------

.. currentmodule:: mrsimulator.spin_system.tensors

.. autoclass:: AntisymmetricTensor
    :show-inheritance:

    .. rubric:: Method Documentation
    .. automethod:: to_freq_dict
    .. automethod:: to_dict_with_units


.. _isotope_api:

Isotope
-------

.. currentmodule:: mrsimulator.spin_system.isotope

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
