
.. _method_documentation:

======
Method
======

Mrsimulator allows users to create custom methods and simulate the NMR spectrum.
At the top level, a :ref:`method_api` object is no different than the pre-built
methods provided within the ``mrsimulator.method.lib`` module.

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
        affine_matrix=[1, 1, 1, 1],
    )

where `name` is an optional method name, `channels` is a list of isotopes used in the
method, `magnetic_flux_density`, `rotor_angle`, and `rotor_frequency` are global
parameters for the method, `spectral_dimension` is the list of SpectralDimension
objects defining the spectral grid, and `affine_matrix` is an optional affine square
matrix.

Although similar to the stock methods from the ``mrsimulator.method.lib`` module, the
above example lacks instructions on how to evaluate frequencies for each spectral dimension.
We pre-defined these instructions for the stock methods for the user's convenience. Here,
we describe how users can write custom instructions.

SpectralDimension
-----------------

A SpectralDimension object is not just a placeholder for defining a spectral grid. It is
also where we define various events---``SpectralEvent`` and ``MixingEvent``, of which the
SpectralEvent is responsible for the NMR frequencies. The syntax for a SpectralDimension
object follows,

.. code-block:: python

    from mrsimulator.method import SpectralEvent, MixingEvent

    SpectralDimension(
        count=512,
        spectral_width=5e4,  # Hz
        reference_offset=10,  # Hz
        origin_offset=4e8,  # Hz
        events=[
            # List of event objects (SpectralEvent and MixingEvent)
            SpectralEvent(name="e0", fraction=0.5),  # fractions are the weights
            # MixingEvent(name="m01"),
            SpectralEvent(name="e1", fraction=0.5),
        ],
    )

where `count`,  `spectral_width`, `reference_offset`, and  `origin_offset` collectively
define the spectral grid, and `events` is a list of spectral and mixing event objects.

The net frequency, :math:`\mathbf{f}_j`, associated with the :math:`j^\text{th}` spectral
dimension is the weighted average of the frequencies from each spectral event within the
dimension,

.. math::
  :label: eq_spectral_average

    \mathbf{f}_j = \sum_{i=0}^{N-1} ~ w_i ~~ \mathbf{e}_i,

where the index :math:`i` spans through the list of spectral events, and :math:`w_i` and
:math:`\mathbf{e}_i` are the weight and corresponding frequency vector from the
:math:`i^\text{th}` spectral event.

In the above example, the average frequency is
:math:`\mathbf{f} = 0.5 \mathbf{e}_0 + 0.5 \mathbf{e}_1`.

.. note::
  Mixing events are not directly involved in spectral frequencies.



Events
------

SpectralEvent
'''''''''''''

A SpectralEvent is where we add instructions on how the frequencies are calculated in mrsimulator.
A generic syntax for the ``SpectralEvent`` follows,

.. code-block:: python

    SpectralEvent(
        fraction=0.5,  # weights w_i
        magnetic_flux_density=4.7,  # T
        rotor_angle=57.735 * 3.1415 / 180,  # rad
        rotor_frequency=10000,  # Hz
        freq_contrib=["Quad2_0", "Quad2_4"],  # frequency contributions list.
        transition_query=[
            {"ch1": {"P": [-3], "D": [0]}},  # A TransitionQuery object
        ],  # transition queries list
    )

Here, `fraction` is the frequency scaling factor for the event and is the same as the weight,
:math:`w_i` in Eq. :eq:`eq_spectral_average`. The attributes `magnetic_flux_density`,
`rotor_angle`, and `rotor_frequency` describe the condition under which frequencies are computed.
These attributes are local to the event, `i.e.`, attributes from a spectral event do not
carry over to the next spectral event. If undefined, the global value from the method attribute
is used for the event.

The attribute `freq_contrib` is a list of frequency contributions allowed during the
event and is used to select specific frequency contributions.
In the above example, the selection only allows the second-order zeroth and fourth-rank
quadrupolar frequency contributions during the event. If undefined, all frequency
contributions are allowed by default. Refer to the :ref:`freq_contrib_api` for the list of
allowed enumerations and :numref:`tb_freq_components` for further details.

The attribute `transition_query` is a list of TransitionQuery objects. These objects query
the SpinSystem objects for a set of allowed spin transitions during the event, `i.e.`, the
ones that satisfy the queries selection criterion. In the above example, we specify a single
TransitionQuery that queries the spin system objects for transitions
that satisfy :math:`p= m_f - m_i = -3` and :math:`d=m_f^2 - m_i^2=0` on channel-1, where
:math:`m_f` and :math:`m_i` are the spin quantum number for the final and initial energy
states involved in a spin-transition. The index `1` in `ch1` is relative to the channels
specified within the method object. In this case, `ch1` refers to ``27Al``.
For details, read the documentation on :ref:`query_doc`.


MixingEvent
'''''''''''
Unlike SpectralEvent, a mixing event is not directly involved in frequency computation. When
a method uses multiple spectral events, each spectral event may query and select a set
of allowed spin transitions. The job of a mixing event is to select which spin
transition from a spectral event, say **e0**, will mix with the spin transitions from the
subsequent spectral event **e1**. As such, mixing events are generally sandwiched between
two spectral events, as follows,

.. code-block:: python

    SpectralDimension(
        events=[
            SpectralEvent(name="e0", fraction=0.5),
            MixingEvent(name="m01", query={"ch1": {"angle": 3.14159, "phase": 0}}),
            SpectralEvent(name="e1", fraction=0.5),
        ],
    )

A MixingEvent object contains the attribute `query`, whose value is a MixingQuery
object. In the above example, the mixing query object queries channel-1, ``27Al``,
for all allowed transitions from spectral events, **e0**, that when rotated by :math:`\pi`
with a phase zero, results in a transition allowed by the spectral event, **e1**. The
resulting pair of transitions form a set of allowed transition pathways.

Examples
--------

**A one-dimension isotropic 3Q-MAS projection**

:math:`\mathbf{\nu}_\text{iso} =  \frac{9}{16}\nu_{3Q} + \frac{7}{16}\nu_{1Q}`

.. code-block:: python

    SpectralDimension(
        events=[
            SpectralEvent(
                fraction=9 / 16, transition_query=[{"ch1": {"P": [-3], "D": [0]}}]
            ),
            SpectralEvent(
                fraction=7 / 16, transition_query=[{"ch1": {"P": [-1], "D": [0]}}]
            ),
        ]
    )

**A one-dimensional Hahn echo**

:math:`\mathbb{p}: +1 \xrightarrow[]{\pi} -1`

.. code-block:: python

    SpectralDimension(
        events=[
            SpectralEvent(fraction=0.5, transition_query=[{"ch1": {"P": [1]}}]),
            MixingEvent(query={"ch1": {"angle": 3.14159, "phase": 0}}),
            SpectralEvent(fraction=0.5, transition_query=[{"ch1": {"P": [-1]}}]),
        ]
    )

**A one-dimensional solid echo**

:math:`\mathbb{p}: -1 \xrightarrow[]{\frac{\pi}{2}} -1`

.. code-block:: python

    SpectralDimension(
        events=[
            SpectralEvent(fraction=0.5, transition_query=[{"ch1": {"P": [-1]}}]),
            MixingEvent(query={"ch1": {"angle": 3.14159 / 2, "phase": 0}}),
            SpectralEvent(fraction=0.5, transition_query=[{"ch1": {"P": [-1]}}]),
        ]
    )

Reference Tables
----------------

.. cssclass:: table-bordered table-striped centered
.. _table_method:
.. list-table:: The attributes of a Method, Method1D, and Method2D object
  :widths: 20 15 65
  :header-rows: 1

  * - Attribute Name
    - Type
    - Description

  * - channels
    - ``List``
    - A *required* list of isotopes given as strings over which the given method applies.
      For example, ``["1H"]``.

  * - magnetic_flux_density
    - ``float``
    - An *optional* float describing the macroscopic magnetic flux density of the applied
      external magnetic field in tesla. For example, ``18.8`` tesla. The default value is
      ``9.4`` tesla.

  * - rotor_frequency
    - ``float``
    - An *optional* float describing the sample rotation frequency in Hz. For example, ``2000`` Hz.
      The default value is ``0`` Hz.

  * - rotor_angle
    - ``float``
    - An *optional* float describing the angle between the sample rotation axis and the external
      magnetic field in radians. The default value is the magic angle,
      ``54.735 * 3.14159 / 180 = 0.955305`` radians.

  * - spectral_dimensions
    - ``List``
    - A list of :ref:`spectral_dim_api` objects describing the spectral dimensions for the method.

  * - affine_matrix
    - ``np.ndarray``
    - A (``n`` x ``n``) affine transformation matrix represented by a numpy array where ``n`` is
      the number of spectral dimensions. If provided, the transformation is applied after running
      a simulation. The default value is ``None`` and no transformation is applied.

  * - simulation
    - CSDM object
    - A CSDM object representing the spectrum simulated by the method. By default, the value is
      ``None``. A value is assigned to this attribute when you run the
      simulation using the :py:meth:`~mrsimulator.Simulator.run` method.

  * - experiment
    - CSDM object
    - An *optional* CSDM object holding an experimental measurement of the method. The default
      value is ``None``


.. cssclass:: table-bordered table-striped centered
.. _table_spectral_dim:
.. list-table:: The attributes of a SpectralDimension object
  :widths: 20 15 65
  :header-rows: 1

  * - Attribute Name
    - Type
    - Description

  * - count
    - ``int``
    - An *optional* integer representing the number of points, :math:`N`, along the spectroscopic
      dimension. For example, ``4096``. The default value is ``1024``.

  * - spectral_width
    - ``float``
    - An *optional* float representing the width, :math:`\Delta x`, of the spectroscopic dimension
      in Hz. For example, ``10e3`` for 10 kHz. The default value is ``25000`` Hz.

  * - reference_offset
    - ``float``
    - An *optional* float representing the reference offset, :math:`x_0`, of the spectroscopic
      dimension in Hz. For example, ``-8000`` Hz. The default value is ``0``.

  * - origin_offset
    - ``float``
    - An optional float representing the origin offset, or Larmor frequency, along the
      spectroscopic dimension in units of Hz. The default value is ``None`` and the origin offset
      is set to the Larmor frequency of isotope from the :attr:`~mrsimulator.Method.channels`
      attribute of the method containing the spectral dimension.

  * - events
    - ``List``
    - An *optional* list of :ref:`event_api` objects used to emulate an experiment.
      The default value is a list with a single **SpectralEvent** with a symmetry_query of
      P=[-1]


.. cssclass:: table-bordered table-striped centered
.. _table_spectral_event:
.. list-table:: The attributes of a SpectralEvent object
  :widths: 20 15 65
  :header-rows: 1

  * - Attribute Name
    - Type
    - Description

  * - magnetic_flux_density
    - ``float``
    - An *optional* float describing the macroscopic magnetic flux density of the applied
      external magnetic field in tesla. For example, ``18.8`` tesla. The default value is
      ``None`` and takes the global magnetic flux density defined by
      :attr:`~mrsimulator.Method.magnetic_flux_density`.

  * - rotor_angle
    - ``float``
    - An *optional* float describing the angle between the sample rotation axis and the external
      magnetic field in radians. The default is ``None`` and takes the global rotor angle defined
      by :attr:`~mrsimulator.Method.rotor_angle`.

  * - rotor_frequency
    - ``float``
    - An *optional* float describing the sample rotation frequency in Hz. For example, ``2000`` Hz.
      The default value is ``None`` and takes the global rotor frequency defined by
      :attr:`~mrsimulator.Method.rotor_frequency`.

  * - freq_contrib
    - ``List``
    - An *optional* list of :ref:`freq_contrib_api` ((object?)) selecting which frequency
      contributions to include when calculating the spectrum. For example,
      ``["Shielding1_0", "Shielding1_2"]``. By default, the list is all frequency enumerations and
      all frequency contributions are calculated.

  * - transition_query
    - ``dict`` or :ref:`transition_api`
    - An *optional* ``dict`` or :ref:`transition_api` selecting transitions active
      during the event. Only these selected transitions will contribute to the net frequency.


.. cssclass:: table-bordered table-striped centered
.. _table_mixing_event:
.. list-table:: The attributes of a MixingEvent object
  :widths: 20 15 65
  :header-rows: 1

  * - Attribute Name
    - Type
    - Description

  * - query
    - ``dict``
    - A mixing_query object selecting a set of transition pathways between two SpectralEvents

..   - The coordinates along each spectral dimension are
..       described with the keywords, *count* (:math:`N`), *spectral_width*
..       (:math:`\nu_\text{sw}`), and *reference_offset* (:math:`\nu_0`). The
..       coordinates are evaluated as,
..
..       .. math
..         \left([0, 1, 2, ... N-1] - \frac{T}{2}\right) \frac{\nu_\text{sw}}{N} + \nu_0
..
..       where :math:`T=N` when :math:`N` is even else :math:`T=N-1`.
