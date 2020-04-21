
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
        >>> site1.isotope.spin
        0.5
        >>> site1.isotope.natural_abundance # in %
        1.11
        >>> site1.isotope.gyromagnetic_ratio # in MHz/T
        10.7084
        >>> site1.isotope.atomic_number
        6
        >>> site1.isotope.quadrupole_moment # in eB
        0.0

    .. rubric:: Method Documentation

    .. automethod:: parse_dict_with_units
    .. automethod:: to_freq_dict
    .. automethod:: to_dict_with_units
