.. _isotopomer_api:

Isotopomer
==========

.. currentmodule:: mrsimulator

.. autoclass:: Isotopomer
    :show-inheritance:

    .. autoattribute:: Zeeman_energy_states
    .. autoattribute:: all_transitions

    .. rubric:: Method Documentation

    .. automethod:: parse_dict_with_units
    .. automethod:: to_freq_dict
    .. automethod:: to_dict_with_units
    .. doctest::

        >>> pprint(isotopomer_1.to_dict_with_units())
        {'abundance': '100 %',
         'sites': [{'isotope': '13C',
                    'isotropic_chemical_shift': '20.0 ppm',
                    'shielding_symmetric': {'eta': 0.5, 'zeta': '10.0 ppm'}}]}

    .. automethod:: get_isotopes
