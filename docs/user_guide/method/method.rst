
.. _method_documentation:

======
Method
======

**mrsimulator**'s organization of the :ref:`spin_sys_api` object and its
component objects, :ref:`site_api`, and :ref:`coupling_api` are easily
understood by anyone familiar with their corresponding physical concepts. The
organization of the :ref:`method_api` object in **mrsimulator** and its related
component objects, however, requires a more detailed explanation of its design.

Overview
--------

An experimental NMR method involves a sequence of rf pulses, free evolution
periods, and sample motion (most commonly magic-angle spinning in solids).
The :ref:`method_api` object in **mrsimulator** only models NMR pulse sequences that
lead to a frequency-domain signal, i.e., a spectrum. The :ref:`method_api`
object is designed to be versatile in its ability to model spectra from
various multi-pulse NMR methods using concepts from the `symmetry pathway
approach <https://doi.org/10.1016/j.pnmrs.2010.11.003>`_ where a pulse
sequence is understood in terms of a set of desired (and undesired) 
*transition pathways*. Each transition pathway is associated with a single
resonance in a multi-dimensional NMR spectrum. The transition pathway signal
encodes information about the spin system interactions in its amplitude and
correlated frequencies. Consider the illustration of a 2D pulse sequence is
shown below, where a desired signal for the method is associated with a
particular transition pathway, :math:`{\hat{A} \rightarrow \hat
{B} \rightarrow \hat{C} \rightarrow \hat{D} \rightarrow \hat
{E} \rightarrow \hat{F}}`.


.. figure:: ../../_static/TransitionPathway.*
    :width: 700
    :alt: figure
    :align: center

Here, the first spectral dimension, i.e., the Fourier transform of the
transition pathway signal as a function of :math:`t_1`, derives its *average 
frequency*, :math:`\overline{\Omega}_1`, from a weighted average of the :math:`\hat{A}`, 
:math:`\hat{B}`, and :math:`\hat{C}` transition frequencies. The second spectral
dimension, i.e., the FT with respect to :math:`t_2`, derives its average frequency, 
:math:`\overline{\Omega}_2`, from a weighted average of the :math:`\hat
{E}`, and :math:`\hat{F}` transition frequencies. Much of the experimental
design and implementation of an NMR method is in identifying the desired
transition pathways and finding ways to acquire their signals while
eliminating all undesired transition pathway signals. 

While NMR measurements take place in the time domain, **mrsimulator** simulates
the corresponding multi-dimensional spectra directly in the frequency domain.
The :ref:`method_api` object in **mrsimulator** needs only
a few details of the NMR pulse sequence to generate the spectrum.  It mimics
the result of the pulse sequence given the desired transition pathways and
their complex amplitudes and average frequencies in each spectroscopic dimension 
of the dataset. To this end, a :ref:`method_api` object is organized according to 
the UML diagram below.  


.. figure:: ../../_static/MethodUML.*
    :width: 700
    :alt: figure
    :align: center


At the heart of a :ref:`method_api` object, assigned to the attribute
``spectral_dimensions``, is an ordered list of :ref:`spectral_dim_api` objects
in the same order as the time evolution dimensions of the experimental NMR
sequence. In each :ref:`spectral_dim_api` object, assigned to the attribute
``events``, is an ordered list of :ref:`event_api` objects, which are divided
into three types: (1) :py:meth:`~mrsimulator.method.SpectralEvent`,
(2) :py:meth:`~mrsimulator.method.ConstantDurationEvent`, and
(3) :py:meth:`~mrsimulator.method.MixingEvent`.  This ordered list
of :ref:`event_api` objects is used to select the desired transition pathways
and determine their average frequency and complex amplitude in the
:py:meth:`~mrsimulator.method.spectral_dimension.SpectralDimension`.  

:py:meth:`~mrsimulator.method.SpectralEvent` and 
:py:meth:`~mrsimulator.method.ConstantDurationEvent` objects are associated 
with excited states of the spin system, with selected transitions evolving 
under the influence of specified Hamiltonian contributions. No coherence 
transfer among transitions or populations can occur in a spectral or 
constant duration event. **mrsimulator** allows the user to select among a 
list of NMR frequency contributions to transitions present during such an 
event in the ``freq_contrib`` attribute holding a list of 
:ref:`enumeration literals<freq_contrib_api>`.  If unspecified, i.e., its value 
is set to ``Null``, a default list holding the enumeration literals for 
all contributions is generated for the event.


.. note::

  All frequency contributions from direct and indirect spin-spin couplings are 
  calculated in the weak-coupling limit in **mrsimulator**.


Additionally, the user can change other measurement attributes during a spectral
or constant duration event: ``rotor_frequency`` or ``rotor_angle``,
``magnetic_flux_density``.  If unspecified, these attributes default to the
values of the identically named global attributes in the :ref:`method_api` object.
Spectral events objects use the ``fraction`` attribute  to calculate the
weighted average frequency for each selected transition pathway during the
spectral dimension.

Inside :py:meth:`~mrsimulator.method.SpectralEvent` and 
:py:meth:`~mrsimulator.method.ConstantDurationEvent` objects, is a
list of :py:meth:`~mrsimulator.method.query.TransitionQuery` objects (*vide infra*) 
which determine which
transitions are "alive" during the event.  :ref:`method_api` objects in
**mrsimulator** are general purpose in the sense that they are designed for
an arbitrary spin system.  That is, a method does not know, in advance, the
energy eigenvalues and eigenstates of the spin system.  Thus, when designing
a :ref:`method_api` object you cannot identify and select a transition through
its initial and final eigenstate quantum numbers.  Instead, transition selection
is done through :py:meth:`~mrsimulator.method.query.TransitionQuery` 
objects during individual spectral or constant duration events.  At
some point, the :ref:`method_api` object uses its 
:py:meth:`~mrsimulator.method.query.TransitionQuery` objects 
to determine the selected transition pathways for a given :ref:`spin_sys_api` 
object as identified by their initial and final eigenstate quantum numbers. 
:py:meth:`~mrsimulator.method.query.TransitionQuery` 
objects hold a list of :py:meth:`~mrsimulator.method.query.SymmetryQuery`
objects which act on specific isotopes present in the, as yet
to be determined, spin system.  A list of isotopes upon which the 
:py:meth:`~mrsimulator.method.query.SymmetryQuery` objects can act are determined 
by the ``channels`` attribute in :ref:`method_api`.  

Inside :py:meth:`~mrsimulator.method.MixingEvent` objects is a 
:py:meth:`~mrsimulator.method.query.MixingQuery` object, which determines the
coherence transfer amplitude between transitions. A 
:py:meth:`~mrsimulator.method.query.MixingQuery` object holds a list of 
:py:meth:`~mrsimulator.method.query.RotationalQuery` objects which act on 
specific isotopes present in the spin system. The list of isotopes upon which 
the :py:meth:`~mrsimulator.method.query.RotationalQuery` objects can act are 
determined by the channels attribute in Method. 

In this guide for designing custom Method objects, we focus first on the
queries objects, i.e., SymmetryQuery and RotationalQuery, and how to use
design them to select the desired transition pathways for a custom method.
The ability to select from a list of frequency contributions by setting
the ``freq_contrib`` attribute in a SpectralEvent of ConstantDuration object 
to a list of :ref:`enumeration literals<freq_contrib_api>` provides another 
means to ensure that your custom Method objects has the right behavior.

Symmetry Query
--------------

Before giving details on how to create a SymmetryQuery object, we need to
review a few key concepts about spin transition symmetry functions. The 
number of possible transition pathways for a spin system depends on the
number of energy eigenstates and the number of spectral and constant duration
events in a method. The number of quantized energy eigenstates for :math:`N`
coupled nuclei is 

.. math::

    \Upsilon_{\left\{ I_1, I_2, \ldots, I_N \right\}} = \prod_{u=1}^N (2 I_u+1),

where :math:`I_u` is the total spin angular momentum of the :math:`u\text
{th}` nucleus and the system of coupled nuclei under consideration is
represented with the notation 
:math:`\left\{ I_1, I_2, \ldots, I_N \right\}`. The transition from quantized
energy level :math:`E_i` to :math:`E_j` is one of 

.. math::

    \mathcal{N}_{\left\{ I_1, I_2, \ldots, I_N \right\}} = \frac{\Upsilon_{\left\{ I_1, I_2, \ldots, I_N \right\}}!}{(\Upsilon_{\left\{ I_1, I_2, \ldots, I_N \right\}}-2)!}

possible transitions between the :math:`\Upsilon` levels.   Here we
count :math:`i  \rightarrow  j` and :math:`j  \rightarrow  i` as different
transitions.  For example, a single spin with angular momentum :math:`I=3/2`,
indicated by :math:`\left\{ I \right\} = \left\{ \tfrac{3}{2} \right\}`, will
have :math:`\Upsilon_{\left\{ 3/2 \right\}} = 2I+1 = 4` energy levels
and :math:`\mathcal{N}_{\left\{ 3/2 \right\}} = 2I(2I+1) = 12` possible NMR
transitions.   A two spin system, :math:`\left\{ I, S \right\} = \left\{ \tfrac
{1}{2}, \tfrac{1}{2} \right\}`, will have 

.. math::

    \Upsilon_{\left\{ 1/2, 1/2 \right\}} = (2I +1) \cdot (2S +1) = 4

 energy levels and

.. math::

    \mathcal{N}_{\left\{ 1/2,1/2 \right\}} =  
    \frac{[(2I +1) \cdot (2S +1)]!}{((2I +1) \cdot (2S +1)-2)!} 
    = \frac{[2 \cdot 2]!}{(2 \cdot 0)!} = 12

possible NMR transitions. For compactness, we will write a transition
(coherence) from state :math:`i` to :math:`j` using the outer product
notation :math:`\ketbra{j}{i}`.  In **mrsimulator**, all simulations are
performed in the high-field limit and further assume that all spin-spin 
couplings are in the weak limit.  In this case, we can identify a transition  
by the quantum numbers of its initial and final Zeeman eigenstate. In the density 
matrix for an ensemble of a given spin system, we could easily identify 
a transition by its row and column indexes.  However, those indexes depend 
on how you have assigned the spins and their eigenstates to those indexes.  
Remember, we need to design the Method object without any details of the 
spin systems upon which they will act.

Selecting Single-Spin Transitions
'''''''''''''''''''''''''''''''''

One way you can select a subset of single-spin transitions, if you don't know the
spin :math:`I` and its associated energy eigenstate quantum numbers, is to request 
all transitions whose spin transition symmetry function, :math:`\text{p}_I` symmetry 
function is :math:`-1`, i.e.,

.. math::
    \text{p}_I(m_f,m_i) = m_f - m_i = -1.

The :math:`\text{p}_I` spin transition symmetry function is also known as the
`coherence order of the transition <https://doi.org/10.1016/0022-2364
(84)90142-2>`_.  

.. note::

    In the high field limit, only single-spin transitions with 
    :math:`{\text{p}_I = \pm 1}` are directly observed.  For a given single-spin
    transition, the signals from :math:`{\text{p}_I = \pm 1}` are complex conjugates 
    of each other, so the convention is to only present the :math:`{\text{p}_I = - 1}`` 
    transition signal in spectra.  

By selecting only single-spin transitions with :math:`\text{p}_I = -1`, you get all
the "observed" transitions from the 
set of all possible transitions.  Similarly, you can use  :math:`\text{p}_I` 
to select any subset of single-spin transitions, such as double-quantum 
(:math:`\text{p}_I = \pm 2`) transitions, 
triple-quantum (:math:`\text{p}_I = \pm 3`) transitions, etc.

Specifying :math:`\text{p}_I` alone, however, is not enough to select
an individual single-spin transition.  However, any individual single-spin 
transition can be identified by a combination of :math:`\text{p}_I` and the 
single-spin transition symmetry function :math:`\text{d}_I`, given by

.. math::

    \text{d}_I(m_i,m_j) =  ~m_j^2 - m_i^2.

Shown below are the values of :math:`\text{p}_I` and :math:`\text{d}_I` 
for all single-spin transitions for :math:`I=1`, :math:`I=3/2` 
and :math:`I=5/2`


.. figure:: ../../_static/SpinOneThreeHalves.*
    :width: 800
    :alt: figure
    :align: center


.. figure:: ../../_static/SpinFiveHalf.*
    :width: 800
    :alt: figure
    :align: center



.. note::

    In the `symmetry pathway approach
    <https://doi.org/10.1016/j.pnmrs.2010.11.003>`_,  the idea of coherence order is extended to form
    a complete set of spin transition symmetry functions, :math:`{\xi}_l
    (i,j)`, given by

    .. math::

        \xi_\ell(i,j) = \bra{j}  \hat{T}_{\ell,0} \ket{j} - \bra{i}  \hat{T}_{\ell,0} \ket{i},

    where the :math:`\hat{T}_{l,0}` are irreducible tensor operators.  The function
    symbol :math:`\xi_\ell(i,j)` is replaced with the lower-case symbols  
    :math:`\mathbb{p}(i,j)`, :math:`\mathbb{d}(i,j)`, :math:`\mathbb{f}
    (i,j)`, :math:`\ldots`, i.e., we follow the spectroscopic sub-shell letter
    designations:

    .. math::

        \begin{array}{cccccccccccccccl}
        \ell = & 0 & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & 10  &11  &12  &13  & \leftarrow \text{numerical value} \\
        \xi_\ell \equiv	& \mathbb{s} &  \mathbb{p} &  \mathbb{d} &  \mathbb{f} &  \mathbb{g} &  \mathbb{h} &  \mathbb{i} & \mathbb{k} &\mathbb{l} & \mathbb{m} & \mathbb{o} & \mathbb{q} & \mathbb{r} &\mathbb{t} & \leftarrow \text{symbol}\\
        \end{array}

    To simplify usage in figures and discussions, we scale the transition symmetry
    functions to integers values according to

    .. math::

        \text{p}(i,j) = \mathbb{p}(i,j), ~~~~~
        \text{d}(i,j) = \sqrt{\frac{2}{3}} \, \mathbb{d}(i,j), ~~~~~
        \text{f}(i,j) = \sqrt{\frac{10}{9}} \, \mathbb{f}(i,j),
        ~~~~~
        \cdots

    The :math:`\ell=0` function is dropped as it always evaluates to zero. For a
    single spin, :math:`I`, a complete set of functions are defined up to 
    :math:`\ell = 2I`. As described in ":ref:`theory`", these functions functions 
    play an important role in evaluating the individual frequency contributions in
    given in :py:meth:`~mrsimulator.method.frequency_contrib.FrequencyEnum` to the
    overall transition frequency. They can also be used to design pulse sequences by
    identifying how different frequency contributions refocus through the
    transition pathways.


For example, for spin :math:`I=1`, the transition :math:`\ketbra{-1}{0}` 
can be selected with :math:`(\text{p}_I,\text{d}_I) = (-1,1)`.  In **mrsimulator**, 
this transition is specified using the SymmetryQuery object, 

.. code-block:: python

    from mrsimulator.method.query import SymmetryQuery

    symm_query = SymmetryQuery(P = [-1], D=[1])


Most notable, :math:`\text{d}_I = 0` for all symmetric transitions, 
:math:`m \rightarrow - m`.  This is particularly useful for quadrupolar 
nuclei, as these transitions are unaffected by the first-order quadrupolar 
coupling frequency contribution.  Thus,  :math:`\ketbra{-\tfrac{1}{2}}{\tfrac{1}{2}}` 
the so-called "central transition" of a quadrupolar nucleus, can be specified 
with the SymmetryQuery object below

.. code-block:: python

    ct_query = SymmetryQuery(P = [-1], D=[0])

Similarly, the symmetric triple quantum transition 
:math:`\ketbra{-\tfrac{3}{2}}{\tfrac{3}{2}}` can be specified using 

.. code-block:: python

    sym_triple_quant_query = SymmetryQuery(P = [-3], D=[0])



Selecting Multi-Spin Transitions
'''''''''''''''''''''''''''''''''


A SymmetryQuery object has only two attributes ``P`` and ``D``, each of which holds a
list of single-spin transition symmetry function values for :math:`

Let's look at an example of creating a SymmetryQuery object that will be associated
with the ``ch1`` attribute of a TransitionQuery object in a SpectralEvent.

.. figure:: ../../_static/CoupledOneHalf.*
    :width: 700
    :alt: figure
    :align: center


Rotational Query
----------------


Mixing events are used to transfer (permute) among transitions and populations,
e.g., :math:`\pi/2` or :math:`\pi` rotations between consecutive spectral or
constant duration events.  For a rotation in a mixing event, the efficiency
associated with the coherence transfer from 

.. math::
    :label: transition

    \ketbra{I, m_f}{I, m_i} \stackrel{\theta_\phi}{\longrightarrow} a(\theta,\phi) \ketbra{I,m_f'}{I,m_i'}

is 

.. math::
    :label: rotation

     a(\theta,\phi) = d_{m_f',m_f}^{(I)}(\theta)d_{m_i',m_i}^{(I)}(\theta)e^
    {-i\Delta p\phi}(i)^{\Delta p}

where :math:`\Delta p = p' - p`.  From this result, we obtain a useful
rule that

.. math::
    :label: piPulseTransition

    \ketbra{m_f}{m_i}  \stackrel{\pi_\phi}{\longrightarrow} \ketbra{-m_f}{-m_i}
    e^{-i\Delta p\phi}(i)^{\Delta p}

The :py:meth:`~mrsimulator.method.MixingEvent` object holds the details of these
rotations in a :py:meth:`~mrsimulator.method.query.MixingQuery` object as a 
:py:meth:`~mrsimulator.method.query.RotationalQuery` object associated with a
``channels`` attribute.

.. code-block:: python

    import numpy as np
    from mrsimulator.method.query import RotationalQuery

    rotation = RotationalQuery(angle = np.pi/2, phase = 0)

It is through :py:meth:`~mrsimulator.method.query.MixingQuery` and 
:py:meth:`~mrsimulator.method.query.TransitionQuery` 
objects that the desired transition pathways are selected and undesired transition 
pathways are eliminated.



SpectralDimension
-----------------

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


:py:meth:`~mrsimulator.method.spectral_dimension.SpectralDimension` has additional 
attributes that have already been
discussed in earlier sections of the documentation.  Notably, ``origin_offset`` 
and ``reference_offset`` are important for converting
the frequency coordinate into a dimensionless frequency ratio coordinate. For
spectra where all the spectral dimensions are associated with single-quantum
transitions on a single isotope, the convention for defining ``origin_offset`` 
and ``reference_offset`` is well established;
the ``origin_offset``, :math:`o_k`, is interpreted as the NMR spectrometer
frequency and  the ``reference_offset`, :math:`b_k`, as the reference
frequency. Given the frequency coordinate, :math:`{X}`, the corresponding
dimensionless-frequency ratio follows,

.. math::
    :label: chemicalShiftDef

    {X}^\text{ratio} = \displaystyle \frac{{X}}{o_k - b_k}.

In the case of multiple quantum dimensions, however, there appear
to be no formal conventions for defining ``origin_offset`` and 
``reference_offset``. 

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

Attribute Summaries
-------------------

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
