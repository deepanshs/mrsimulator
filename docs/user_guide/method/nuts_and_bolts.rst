
===========================
Method (For advanced users)
===========================

Mrsimulator allows users to create custom methods and simulate the NMR spectrum.
At the top level, a Method object is no different than the stock methods provided
within the mrsimulator methods library.

A generic setup for a custom method (similar to the stock method) follows,

.. code-block:: python

    from mrsimulator.method import Method, SpectralDimension

    my_method = Method(
        name="my_method",
        channels=["27Al", "13C"],  # list of isotopes
        magnetic_flux_density=4.7,  # T
        rotor_angle=57.735 * 3.1415 / 180,  # rad
        rotor_frequency=10000,  # Hz
        spectral_dimensions=[
            SpectralDimension(count=512, spectral_width=50000),  # dimension-0
            SpectralDimension(count=256, spectral_width=10000),  # dimension-1
        ],
    )

where `name` is an optional method name, `channels` are a list of isotopes used in the
method, `magnetic_flux_density`, `rotor_angle`, and `rotor_frequency` are the global
parameters, and `spectral_dimension` is the list of SpectralDimension objects defining
the spectral grid.
Although similar to the stock method from the ``mrsimulator.methods`` library, the
above example lacks instructions on how to evaluate frequencies for each spectral dimension.
We pre-defined these instructions for the stock methods for the user's convenience. Here,
we describe how users can create custom methods.

SpectralDimension
-----------------

A SpectralDimension object is not just a placeholder for defining a spectral grid. It is
also where we define various events, of which, SpectralEvent is responsible for the
NMR frequencies. The average frequency for a given spectral dimension is

.. math::
    \mathbf{f}_j = \sum_{i=0}^{N-1} ~ w_i ~~ \mathbf{e}_i,

where the index :math:`i` spans through the list of spectral events, :math:`w_i` and
:math:`\mathbf{e}_i` are the weight and corresponding frequency vector from the
:math:`i^\text{th}` spectral event, and :math:`\mathbf{f}_j` is the vector of average
frequencies along the :math:`j^\text{th}` spectral dimension.

A syntax for a SpectralDimension object follows,

.. code-block:: python

    from mrsimulator.method import SpectralEvent

    SpectralDimension(
        count=512,
        spectral_width=500,
        reference_offset=10,
        events=[
            SpectralEvent(name="A", fraction=0.5),  # fractions are the weights
            SpectralEvent(name="B", fraction=0.5),  # A list of events
        ],
    )

Here, the average frequency along the spectral dimension is
:math:`\mathbf{f} = 0.5 \mathbf{e}_0 + 0.5 \mathbf{e}_1`.


Events
------

SpectralEvent
'''''''''''''

The SpectralEvent attributes instruct the objects on how and under what conditions
to compute the NMR frequencies, :math:`\mathbf{e}_i`.

The syntax for ``SpectralEvent`` follows,

.. code-block:: python

    SpectralEvent(
        fraction=0.5,
        magnetic_flux_density=4.7,  # T
        rotor_angle=57.735 * 3.1415 / 180,  # rad
        rotor_frequency=10000,
        freq_contrib=["Quad2_0", "Quad2_4"],  # frequency contributions list.
        transition_query=[{"ch1": {"P": [-3], "D": [0]}}],  # transition queries list
    )

Here, `fraction`, :math:`w_i` is the frequency scaling factor for the event.
The attributes `magnetic_flux_density`, `rotor_angle`, and `rotor_frequency` define the
condition under which frequencies are computed. These attributes are local to the event.
If undefined, the global method parameter is used instead.

The attribute `freq_contrib` is a list of frequency contributions allowed during the
event and is used to toggle frequency contributions.
In the above example, the selection only allows the second-order zeroth and fourth-rank
quadrupolar frequency contributions during the event. If undefined, all frequency
contributions are allowed by default. Refer to the :ref:`freq_contrib_api` for the list
of frequency contributions.

The attribute `transition_query` is a list of transition query objects which include
instructions on how to query and select spin transition(s) during the event. The above
example instructs the event to query the spin system objects for transitions that
satisfy :math:`p= m_f - m_i = -3` and :math:`d=m_f^2 - m_i^2=0` for channel-1, where
:math:`m_f` and :math:`m_i` are the spin quantum number for the final and initial energy
states involved in a spin-transition. The index `1` in `ch1` is relative to the channels
specified within the method objects. In this case, `ch1` refers to ``27Al``.
Read more on the transition query.


MixingEvent
'''''''''''
Unlike SpectralEvent, a mixing event is not involved in frequency computation. When
a method uses multiple spectral events, each spectral event may query and select a set
of allowed spin transitions. The job of a mixing event is to select which spin
transition from a spectral event, say **A**, will mix with the spin transitions from
spectral event **B**. As such, mixing events are generally sandwiched between two spectral
events, as follows,

.. code-block:: python

    from mrsimulator.method import MixingEvent

    SpectralDimension(
        events=[
            SpectralEvent(name="A", fraction=0.5),
            MixingEvent(mixing_query={"ch1": {"tip_angle": 3.14159, "phase": 0}}),
            SpectralEvent(name="B", fraction=0.5),
        ],
    )

The MixingEvent object contains the attribute `mixing_query`, whose value is a mixing
query object. In the above example, the mixing query object queries channel-1, ``27Al``,
for all allowed transitions from the previous spectral events, **A**, that when rotated
by :math:`\pi` with a phase zero, results in a transition allowed by the subsequent
spectral event, **B**. The resulting pair of transitions form a set of allowed transition
pathways.

Examples
--------

**A one-dimension isotropic 3Q-MAS projection**

:math:`\mathbf{\nu}_\text{iso} =  \frac{9}{16}\nu_{3Q} + \frac{7}{16}\nu_{1Q}`

.. code-block:: python

    SpectralDimension(
        events=[
            SpectralEvent(fraction=9 / 16, transition_query=[{"ch1": {"P": [-3], "D": [0]}}]),
            SpectralEvent(fraction=7 / 16, transition_query=[{"ch1": {"P": [-1], "D": [0]}}]),
        ]
    )

**A one-dimensional Hahn echo**

:math:`\mathbb{p}: +1 \xrightarrow[]{\pi} -1`

.. code-block:: python

    SpectralDimension(
        events=[
            SpectralEvent(fraction=0.5, transition_query=[{"ch1": {"P": [1]}}]),
            MixingEvent(mixing_query={"ch1": {"tip_angle": 3.14159, "phase": 0}}),
            SpectralEvent(fraction=0.5, transition_query=[{"ch1": {"P": [-1]}}]),
        ]
    )

**A one-dimensional solid echo**

:math:`\mathbb{p}: -1 \xrightarrow[]{\frac{\pi}{2}} -1`

.. code-block:: python

    SpectralDimension(
        events=[
            SpectralEvent(fraction=0.5, transition_query=[{"ch1": {"P": [-1]}}]),
            MixingEvent(mixing_query={"ch1": {"tip_angle": 3.14159 / 2, "phase": 0}}),
            SpectralEvent(fraction=0.5, transition_query=[{"ch1": {"P": [-1]}}]),
        ]
    )
