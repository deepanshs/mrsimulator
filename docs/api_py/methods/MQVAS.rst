
Multi-quantum variable angle spinning method
--------------------------------------------

.. currentmodule:: mrsimulator.methods

.. autofunction:: MQVAS

**Example of MQ-MAS method for I=3/2**

For I=3/2, the p coherence pathway follows :math:`-3 \rightarrow -1`, where the
resonance from the first and the second block arises from p=-3 and p=-1, respectively.
Use the `transition_query` attribute to specify the `p` values for the first and
second event block. For central-transition selective transition, set d=0.

In the following example, the `magnetic_flux_density` and `rotor_angle` are set as
global parameters. The `transition_query` is defined for each event, where the first
event (indirect dimension) is p=-3, and the second event (direct dimension) is p=-1.

.. doctest::

    >>> from mrsimulator.methods import MQVAS
    >>> Bloch_CT_method = MQVAS(
    ...     channels=['23Na'],
    ...     magnetic_flux_density=9.4, # in T
    ...     rotor_angle=0.95531, # in rad (MAS)
    ...     spectral_dimensions=[
    ...         # The first spectral dimension block is the indirect-dimension
    ...         {
    ...             "count": 512,
    ...             "spectral_width": 5e4, # in Hz
    ...             "reference_offset": -8e3, # in Hz
    ...             "events": [{
    ...                 "transition_query": {"P": [-3], "D": [0]}
    ...             }]
    ...         },
    ...         # The second spectral dimension block is the direct-dimension
    ...         {
    ...             "count": 1024,
    ...             "spectral_width": 2e4, # in Hz
    ...             "reference_offset": -3e3, # in Hz
    ...             "events": [{
    ...                 "transition_query": {"P": [-1], "D": [0]}
    ...             }]
    ...         }
    ...     ]
    ... )


**Example of COASTER method**

A correlation of anisotropies separated through echo refocusing (COASTER) method is the
same as the MQMAS, except the rotor angle is :math:`70.12^\circ` instead of the
magic-angle, and the p=3 for the first the block.

.. doctest::

    >>> from mrsimulator.methods import MQVAS
    >>> Bloch_CT_method = MQVAS(
    ...     channels=['23Na'],
    ...     magnetic_flux_density=9.4, # in T
    ...     rotor_angle=1.223824, # in rad
    ...     spectral_dimensions=[
    ...         # The first spectral dimension block is the indirect-dimension
    ...         {
    ...             "count": 512,
    ...             "spectral_width": 5e4, # in Hz
    ...             "reference_offset": -8e3, # in Hz
    ...             "events": [{
    ...                 "transition_query": {"P": [3], "D": [0]}
    ...             }]
    ...         },
    ...         # The second spectral dimension block is the direct-dimension
    ...         {
    ...             "count": 1024,
    ...             "spectral_width": 2e4, # in Hz
    ...             "reference_offset": -3e3, # in Hz
    ...             "events": [{
    ...                 "transition_query": {"P": [-1], "D": [0]}
    ...             }]
    ...         }
    ...     ]
    ... )

**Example of Switched-angle spinning (SAS) method**

A Switched-angle spinning (SAS) method is a correlation of p=-1 frequencies at two
different angles. The following example shows the SAS method, where the two angles
are :math:`70.12^\circ` and :math:`54.74^\circ`, respectively.

.. doctest::

    >>> from mrsimulator.methods import MQVAS
    >>> Bloch_CT_method = MQVAS(
    ...     channels=['23Na'],
    ...     magnetic_flux_density=9.4, # in T
    ...     spectral_dimensions=[
    ...         # The first spectral dimension block is the indirect-dimension
    ...         {
    ...             "count": 512,
    ...             "spectral_width": 5e4, # in Hz
    ...             "reference_offset": -8e3, # in Hz
    ...             "events": [{
    ...                 "rotor_angle": 1.223824, # 70.12 deg in rads
    ...                 "transition_query": {"P": [-1], "D": [0]}
    ...             }]
    ...         },
    ...         # The second spectral dimension block is the direct-dimension
    ...         {
    ...             "count": 1024,
    ...             "spectral_width": 2e4, # in Hz
    ...             "reference_offset": -3e3, # in Hz
    ...             "events": [{
    ...                 "rotor_angle": 0.95531, # magic-angle in rad
    ...                 "transition_query": {"P": [-1], "D": [0]}
    ...             }]
    ...         }
    ...     ]
    ... )
