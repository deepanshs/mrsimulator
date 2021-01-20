.. _spin_sys_api:

SpinSystem
==========

.. currentmodule:: mrsimulator

.. autoclass:: SpinSystem
    :show-inheritance:

    .. rubric:: Method Documentation

    .. automethod:: get_isotopes
    .. automethod:: zeeman_energy_states
    .. automethod:: all_transitions
    .. automethod:: parse_dict_with_units
    .. automethod:: json
    .. doctest::

        >>> pprint(spin_system_1.json())
        {'abundance': '100 %',
         'sites': [{'isotope': '13C',
                    'isotropic_chemical_shift': '20.0 ppm',
                    'shielding_symmetric': {'eta': 0.5, 'zeta': '10.0 ppm'}}]}
