.. _coupling_api:

Coupling
========

.. currentmodule:: mrsimulator

.. autoclass:: Coupling
    :show-inheritance:

    .. rubric:: Method Documentation

    .. automethod:: parse_dict_with_units
    .. automethod:: json
    .. doctest::

        >>> pprint(coupling1.json())
        {'isotropic_j': '20.0 Hz',
         'j_symmetric': {'eta': 0.5, 'zeta': '10.0 Hz'},
         'site_index': [0, 1]}
