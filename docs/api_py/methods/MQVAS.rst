
Multi-quantum variable angle spinning method
--------------------------------------------

.. currentmodule:: mrsimulator.methods

.. autofunction:: Method2D

**Example of MQ-MAS method for I=3/2**

In an MQ-MAS experiment for spin I=3/2, the p coherence pathway follows
:math:`-3 \rightarrow -1`, where the resonance from the first and the second dimensions
arises from p=-3 and p=-1, respectively. Use the `transition_query` attribute to
specify the `p` values for the events within the two dimensions. For central-transition
selective transition, set d=0.

In the following example, the `magnetic_flux_density` and `rotor_angle` are set as
global parameters. The `transition_query` is defined for each event, where the first
event (indirect dimension) is p=-3, and the second event (direct dimension) is p=-1.

.. doctest::

    >>> from mrsimulator.methods import Method2D
    >>> mq_mas = Method2D(
    ...     channels=['23Na'],
    ...     magnetic_flux_density=9.4, # in T
    ...     rotor_angle=0.95531, # in rad (magic-angle)
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

**Example of ST-MAS method for I=5/2**

Satellite transition correlation method is a 2D method for corelating the satellite
transitions with the central transition. You may create a `stmas` method from ``Method2D``
class similar to the `mq_mas` method in the previous example. The main
difference is setting up the transition query. Instead of a querying for the symmetric
triple-quantum transition, in the case of mqmas method, we will query for satellite
transitions.

For selecting central transition, use p=-1, d=0.

To select the inner staellites, use p=1 and d=[-2, 2]. The given query will select all
transitions with :math:`m_f - m_i=1` and :math:`m_f^2 - m_i^2 = \pm2`, which corresponds
to the inner satellite transitions,
:math:`\left|-\frac{3}{2}\right\rangle \rightarrow \left|-\frac{1}{2}\right\rangle` and
:math:`\left|\frac{1}{2}\right\rangle \rightarrow \left|\frac{3}{2}\right\rangle`

Similarly, to select the outter staellites, use p=-1 and d=[-4, 4]. The given query will
select all transitions with :math:`m_f - m_i=-1` and :math:`m_f^2 - m_i^2 = \pm4`, which
corresponds to,
:math:`\left|\frac{5}{2}\right\rangle \rightarrow \left|\frac{3}{2}\right\rangle` and
:math:`\left|-\frac{3}{2}\right\rangle \rightarrow \left|-\frac{5}{2}\right\rangle`

You may add any or all CT :math:`\rightarrow` CT, ST1 :math:`\rightarrow` CT, or
ST2 :math:`\rightarrow` CT query to the transition query selector. In following example,
we add CT :math:`\rightarrow` CT and ST1 :math:`\rightarrow` CT transition query.


.. doctest::

    >>> from mrsimulator.methods import Method2D
    >>> st_mas = Method2D(
    ...     channels=['27Al'],
    ...     magnetic_flux_density=9.4, # in T
    ...     rotor_angle=1.223824, # in rad (magic-angle)
    ...     spectral_dimensions=[
    ...         # The first spectral dimension block is the indirect-dimension
    ...         {
    ...             "count": 512,
    ...             "spectral_width": 5e4, # in Hz
    ...             "reference_offset": -8e3, # in Hz
    ...             "events": [{
    ...                 "transition_query": {"P": [-1], "D": [0]}
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
magic-angle, and the p=3 for the first the dimension.

.. doctest::

    >>> from mrsimulator.methods import Method2D
    >>> coaster = Method2D(
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

    >>> from mrsimulator.methods import Method2D
    >>> sas = Method2D(
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
