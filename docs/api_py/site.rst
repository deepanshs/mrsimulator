.. _site_api:

Site
====

.. currentmodule:: mrsimulator

.. autoclass:: Site
    :show-inheritance:

    .. rubric:: Method Documentation

    .. automethod:: parse_dict_with_units
    .. automethod:: json
    .. doctest::

        >>> pprint(site1.json())
        {'isotope': '13C',
         'isotropic_chemical_shift': '20.0 ppm',
         'shielding_symmetric': {'eta': 0.5, 'zeta': '10.0 ppm'}}
