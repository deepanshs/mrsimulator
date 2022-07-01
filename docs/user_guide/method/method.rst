
.. _method_documentation:

======
Method
======

While **mrsimulator**'s organization of the :ref:`spin_sys_api` object and its
composite objects, :ref:`site_api`, and :ref:`coupling_api` are easily
understood by anyone familiar with the underlying physical concepts, the
organization of the :ref:`method_api` object in **mrsimulator** and its related
composite objects requires a more detailed explanation of its design.  This
section assumes that you are already familar with topics covered in the
introduction sections :ref:`getting_started`, 
:ref:`introduction_isotopomers_example`, and :ref:`fitting_example`.

Overview
--------

An experimental NMR method involves a sequence of rf pulses, free evolution
periods, and sample motion. The :ref:`method_api` object in **mrsimulator** 
models the spectrum from an NMR pulse sequences. The :ref:`method_api`
object is designed to be versatile in its ability to model spectra from various
multi-pulse NMR methods using concepts from the `symmetry pathway approach
<https://doi.org/10.1016/j.pnmrs.2010.11.003>`_ where a pulse sequence is
understood in terms of a set of desired (and undesired) 
*transition pathways*. Each transition pathway is associated with a single
resonance in a multi-dimensional NMR spectrum. The transition pathway signal
encodes information about the spin system interactions in its amplitude and
correlated frequencies. Consider the illustration of a 2D pulse sequence
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
frequency*, :math:`\overline{\Omega}_1`, from a weighted average of
the :math:`\hat{A}`, :math:`\hat{B}`, and :math:`\hat{C}` transition
frequencies. The second spectral dimension, i.e., the FT with respect
to :math:`t_2`, derives its average frequency, :math:`\overline{\Omega}_2`, 
from a weighted average of the :math:`\hat{E}`, and :math:`\hat{F}` transition 
frequencies. Much of the experimental design and implementationof an NMR method 
is in identifying the desired transition pathways and finding
ways to acquire their signals while eliminating all undesired transition
pathway signals. 

While NMR measurements take place in the time domain, **mrsimulator** simulates
the corresponding multi-dimensional spectra directly in the frequency domain.
The :ref:`method_api` object in **mrsimulator** needs only a few details of the
NMR pulse sequence to generate the spectrum.  It mimics the result of the pulse
sequence given the desired transition pathways and their complex amplitudes and
average frequencies in each spectroscopic dimension of the dataset. To this
end, a :ref:`method_api` object is organized according to the UML diagram
below.  


.. figure:: ../../_static/MethodUML.*
    :width: 700
    :alt: figure
    :align: center

.. note::

  In UML (Unified Modeling Language) diagrams, each class is represented with 
  a box that contains two compartments.  The top compartment contains the name 
  of the class, and the bottom compartment contains the attributes of the class.
  Default attribute values are shown as assignments. A composition 
  is depicted as a binary association decorated with a filled black diamond. 
  Inheritance is shown as a line with a hollow triangle as an arrowhead.


At the heart of a :ref:`method_api` object, assigned to its attribute
``spectral_dimensions``, is an ordered list of :ref:`spectral_dim_api` objects
in the same order as the time evolution dimensions of the experimental NMR
sequence. In each :ref:`spectral_dim_api` object, assigned to the attribute
``events``, is an ordered list of :ref:`event_api` objects, which are divided
into three types: (1) :py:meth:`~mrsimulator.method.SpectralEvent`,
(2) :py:meth:`~mrsimulator.method.DelayEvent`, and
(3) :py:meth:`~mrsimulator.method.MixingEvent`.  This ordered list
of :ref:`event_api` objects is used to select the desired transition pathways
and determine their average frequency and complex amplitude in the
:py:meth:`~mrsimulator.method.spectral_dimension.SpectralDimension`.  

.. warning::

  DelayEvent objects are not available in version 0.7 of **mrsimulator**.


:py:meth:`~mrsimulator.method.SpectralEvent` and 
:py:meth:`~mrsimulator.method.DelayEvent` objects define which 
transitions are "alive" during the event and under which
transition-dependent frequency contributions they evolve. No coherence 
transfer among transitions or populations occurs in a spectral or delay event.
The transition-dependent frequency contributions during an Event are
selected from a list of  :ref:`enumeration literals<freq_contrib_api>` 
and placed in the ``freq_contrib`` attribute of the Event.  If unspecified, 
i.e., the value of ``freq_contrib`` is set to ``Null``, a default list holding 
the enumeration literals for *all* contributions is generated for the event.


.. note::

  All frequency contributions from direct and indirect spin-spin couplings are 
  calculated in the weak-coupling limit in **mrsimulator**.


Additionally, the user can affect transition frequencies during a spectral
or delay event by changing other measurement attributes : ``rotor_frequency`` 
or ``rotor_angle``, ``magnetic_flux_density``.  If unspecified, these attributes 
default to the values of the identically named global attributes in the 
:ref:`method_api`
object. SpectralEvents objects use the ``fraction`` attribute  to calculate
the weighted average frequency for each selected transition pathway during the
spectral dimension.

Inside :py:meth:`~mrsimulator.method.SpectralEvent` and 
:py:meth:`~mrsimulator.method.DelayEvent` objects, is a list
of :py:meth:`~mrsimulator.method.query.TransitionQuery` objects (*vide infra*)
which determine which transitions are "alive" during the
event. :ref:`method_api` objects in
**mrsimulator** are general purpose in the sense that they are designed for an
arbitrary spin system.  That is, a method does not know, in advance, the
energy eigenvalues and eigenstates of the spin system.  Thus, when designing
a :ref:`method_api` object you cannot identify and select a transition
through its initial and final eigenstate quantum numbers.  Instead,
transition selection is done
through :py:meth:`~mrsimulator.method.query.TransitionQuery` objects during
individual spectral or delay events.  It is only during a simulation that
the :ref:`method_api` object uses its 
:py:meth:`~mrsimulator.method.query.TransitionQuery` objects to determine the
selected transition pathways for a given :ref:`spin_sys_api` object by 
the initial and final eigenstate quantum numbers of each transition. 
:py:meth:`~mrsimulator.method.query.TransitionQuery` objects hold a list
of :py:meth:`~mrsimulator.method.query.SymmetryQuery` objects which act on
specific isotopes in the, as yet to be determined, spin system.  A list of
specific isotopes upon which the 
:py:meth:`~mrsimulator.method.query.SymmetryQuery` objects act are determined by
the ``channels`` attribute in :ref:`method_api`.  

Inside :py:meth:`~mrsimulator.method.MixingEvent` objects is a 
:py:meth:`~mrsimulator.method.query.MixingQuery` object, which determines the
coherence transfer amplitude between transitions. A 
:py:meth:`~mrsimulator.method.query.MixingQuery` object holds a list of 
:py:meth:`~mrsimulator.method.query.RotationalQuery` objects which act on
specific isotopes present in the spin system. As before, the list of isotopes
upon which the :py:meth:`~mrsimulator.method.query.RotationalQuery` objects
act are determined by the ``channels`` attribute in Method. 

In this guide to designing custom Method objects, we focus first on the queries
objects, i.e., TransitionQuery and MixingQuery, and how to use them to select
the desired transition pathways for a custom method. Then we examine how
transitions frequencies in the desired transition pathways can be selected from
a list of frequency contributions using the ``freq_contrib`` attribute of a
SpectralEvent of DelayEvent object. The ability to select 
:ref:`frequency contributions<freq_contrib_api>` can often reduce the number of
events needed in the design of your custom Method object.

Theoretical Background
----------------------

Before giving details on how to create a custom Method object, we need 
to review a few key concepts about spin transitions and 
*transition symmetry functions*

The number of quantized energy eigenstates for :math:`N` coupled nuclei is 

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

possible NMR transitions. We write a transition (coherence) from 
state :math:`i` to :math:`j` using the outer product
notation :math:`\ketbra{j}{i}`.  In **mrsimulator**, all simulations are
performed in the high-field limit and further assume that all spin-spin 
couplings are in the weak limit.  

To write a custom Method in *mrsimulator*, you need to determine the
desired transition pathways, and then select the desired transitions during
each SpectralEvent or DelayEvent.  Keep in mind, however, that Method 
objects are designed without any details of the spin systems upon which 
they will act.  For example, in the density matrix of a spin system ensemble,
one could easily identify a transition by its row and column indexes.  
However, those indexes depend on the spin system and
how the spins and their eigenstates have been assigned to those indexes.
Therefore, we can not use such another approach for selecting transitions.

Spin Transition Symmetry Functions
''''''''''''''''''''''''''''''''''

One way you can select a subset of single-spin transitions, even if you don't
know the spin quantum number :math:`I` and its associated energy eigenstate
quantum numbers, is to request all transitions whose single-spin transition
symmetry function, :math:`\text{p}_I` symmetry function is :math:`-1`, i.e.,

.. math::
    \text{p}_I(m_f,m_i) = m_f - m_i = -1.

The :math:`\text{p}_I` single-spin transition symmetry function is also known as
the single-spin `coherence order of the transition
<https://doi.org/10.1016/0022-2364(84)90142-2>`_.  

.. note::

    In the high field limit, only single-spin transitions with 
    :math:`{\text{p}_I = \pm 1}` are directly observed.  Since the
    evolution frequencies of the :math:`\ketbra{j}{i}` and 
    :math:`\ketbra{i}{j}` transitions are equal in magnitude but opposite 
    in sign, the convention is to only present the :math:`{\text{p}_I = - 1}`` 
    transition resonances in single-quantum spectra.  

By selecting only single-spin transitions with :math:`\text{p}_I = -1`, you get
all the "observed" transitions from the set of all possible transitions.
Similarly, you can use  :math:`\text{p}_I` to select any subset of single-spin
transitions, such as double-quantum :math:`(\text{p}_I = \pm 2)` transitions,
triple-quantum :math:`(\text{p}_I = \pm 3)` transitions, etc.

Specifying :math:`\text{p}_I` alone is not enough to select an individual
single-spin transition.  However, any individual single-spin transition can be
identified by a combination of :math:`\text{p}_I` and the single-spin
transition symmetry function :math:`\text{d}_I`, given by

.. math::

    \text{d}_I(m_i,m_j) =  ~m_j^2 - m_i^2.

You can verify this from the values of :math:`\text
{p}_I` and :math:`\text{d}_I` for all single-spin transitions
for :math:`I=1`, :math:`I=3/2` and :math:`I=5/2` shown below.  Note
that :math:`\text{d}_I = 0` for all transitions in a :math:`I=1/2` nucleus.


.. figure:: ../../_static/SpinOneThreeHalves.*
    :width: 600
    :alt: figure
    :align: center


.. figure:: ../../_static/SpinFiveHalf.*
    :width: 650
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
    :math:`\ell = 2I`.
    

    For weakly coupled nuclei, we define the transition symmetry functions

    .. math::

      \xi_{\ell_1,\ell_2, \ldots, \ell_n} (i,j) = 
      \left \langle j \right|\hat{T}_{\ell_1,0}({\bf I}_1)\hat{T}_{\ell_2,0}({\bf I}_2)\ldots\hat{T}_{\ell_n,0}({\bf I}_n) \left|j \right \rangle
      - 
      \left \langle i \right|\hat{T}_{\ell_1,0}({\bf I}_1)\hat{T}_{\ell_2,0}({\bf I}_2)\ldots\hat{T}_{\ell_n,0}({\bf I}_n) \left|i \right \rangle

    Replacing the symmetry function symbol using sub-shell letter designations becomes 
    more cumbersome in this case.  When the :math:`\ell` are zero on all nuclei except one,  
    we identify these functions as

    .. math::

      \begin{array}{cccc}
      \mathbb{p}_1 = \xi_{1,0, \ldots, 0} (i,j), &
      \mathbb{p}_2 = \xi_{0,1, \ldots, 0} (i,j), &
      \ldots, &
      \mathbb{p}_n = \xi_{0,0, \ldots, 1} (i,j),\\
      \\
      \mathbb{d}_1 = \xi_{2, 0, \ldots, 0} (i,j), &
      \mathbb{d}_2 = \xi_{0,2, \ldots, 0} (i,j), &
      \ldots, &
      \mathbb{d}_n = \xi_{0,0, \ldots, 2} (i,j), \\
      \\
      \mathbb{f}_1 = \xi_{3, 0, \ldots, 0} (i,j), &
      \mathbb{f}_2 = \xi_{0,3, \ldots, 0} (i,j), &
      \ldots, &
      \mathbb{f}_n = \xi_{0,0, \ldots, 3} (i,j), \\
      \vdots & \vdots &  & \vdots
      \end{array}

    For weakly coupled homonuclear spins it is also convenient to define 

    .. math::

      \begin{array}{c}
      \mathbb{p}_{1,2,\ldots,n} =  \mathbb{p}_{1} 
      + \mathbb{p}_{2} + \cdots \mathbb{p}_{n} \\
      \\
      \mathbb{d}_{1,2,\ldots,n} =  \mathbb{d}_{1} 
      + \mathbb{d}_{2} + \cdots \mathbb{d}_{n} \\
      \\
      \mathbb{f}_{1,2,\ldots,n} =  \mathbb{f}_{1} 
      + \mathbb{f}_{2} + \cdots \mathbb{f}_{n} \\
      \vdots
      \end{array}


    When the :math:`\ell` are zero on all nuclei except two, then we identify 
    these functions using a combination of sub-shell letter designations, e.g.,

    .. math::

      \begin{array}{cccc}
      (\mathbb{pp})_{1,2} = \xi_{1,1,0, \ldots, 0} (i,j), &
      (\mathbb{pp})_{1,3} = \xi_{1,0,1, \ldots, 0} (i,j), &
      \ldots, &
      (\mathbb{pp})_{1,n} = \xi_{1,0,0, \ldots, 1} (i,j),\\
      \\
      (\mathbb{pd})_{1,2} = \xi_{1, 2, 0, \ldots, 0} (i,j), &
      (\mathbb{pd})_{1,3} = \xi_{1,0,2 \ldots, 0} (i,j), &
      \ldots, &
      (\mathbb{pd})_{1,n} = \xi_{1,0, \ldots, 2} (i,j), \\
      \\
      (\mathbb{dp})_{1,2} = \xi_{2, 1, 0, \ldots, 0} (i,j), &
      (\mathbb{dp})_{1,3} = \xi_{2 ,0, 1 \ldots, 0} (i,j), &
      \ldots, &
      (\mathbb{dp})_{1,n} = \xi_{2, 0, \ldots, 1} (i,j), \\
      \vdots & \vdots &  & \vdots
      \end{array}

    As described in ":ref:`theory`", these transition symmetry functions 
    play an important role in evaluating the individual frequency 
    contributions in given in 
    :py:meth:`~mrsimulator.method.frequency_contrib.FrequencyEnum` to the
    overall transition frequency. They also aid in pulse sequence 
    design by identifying how different frequency contributions 
    refocus through the transition pathways.



Single-Spin Transition Queries
------------------------------

Based on the review above, we now know that the spin :math:`I=1`, the 
transition :math:`\ketbra{-1}{0}` can be selected with 
:math:`(\text{p}_I,\text{d}_I) = (-1,1)`.  In **mrsimulator**, 
this transition is selected during a SpectralEvent using the SymmetryQuery 
and TransitionQuery objects, 

.. plot::
    :context: reset

    from mrsimulator.method.query import SymmetryQuery, TransitionQuery
    from mrsimulator.method import SpectralEvent

    symm_query = SymmetryQuery(P=[-1], D=[1])
    trans_query = TransitionQuery(ch1=symm_query)
    event = SpectralEvent(fraction=1, transition_query=[trans_query])

In the example above, the SymmetryQuery instance is created and assigned to the
``ch1`` attribute of a TransitionQuery, so that it acts on the first isotope in
the list assigned to the ``channels`` attribute  of the Method object.  This
TransitionQuery instance is then added to a list assigned to the 
``transition_queries`` attribute of a SpectralEvent which can be added 
to an ordered list of Events in the ``events`` attribute of a SpectralDimension 
object, as shown in the code below.

.. plot::
    :context: close-figs

    from mrsimulator import Site, Coupling, SpinSystem, Simulator
    from mrsimulator import Method, SpectralDimension
    from mrsimulator import signal_processor as sp
    import matplotlib.pyplot as plt
    import numpy as np

    deuterium = Site(
      isotope="2H",
      isotropic_chemical_shift=10,  # in ppm
      shielding_symmetric={"zeta":-80, "eta":0.25},  # zeta in ppm
      quadrupolar={"Cq":10e3, "eta":0.0,"alpha":0, "beta":np.pi/2, "gamma":0}
    )
    spin_system = SpinSystem(sites=[deuterium])

    method_both_transitions = Method(
      channels=["2H"],
      magnetic_flux_density=9.4,  # in T
      spectral_dimensions=[
        SpectralDimension(
          count=512,
          spectral_width=40000,  # in Hz
          events=[SpectralEvent(fraction=1, 
            transition_query=[{"ch1":{"P":[-1]}}])])
          ]
      )

    method_transition1 = Method(
      channels=["2H"],
      magnetic_flux_density=9.4,  # in T
      spectral_dimensions=[
        SpectralDimension(
          count=512,
          spectral_width=40000,  # in Hz
          events=[SpectralEvent(fraction=1, 
            transition_query=[{"ch1":{"P":[-1], "D":[1]}}])])
          ]
      )

    method_transition2 = Method(
      channels=["2H"],
      magnetic_flux_density=9.4,  # in T
      spectral_dimensions=[
        SpectralDimension(
          count=512,
          spectral_width=40000,  # in Hz
          events=[SpectralEvent(fraction=1, 
            transition_query=[{"ch1":{"P":[-1], "D":[-1]}}])])
          ]
      )

    sim = Simulator(spin_systems = [spin_system],
      methods=[method_both_transitions,method_transition1,method_transition2])
    sim.run()

    processor = sp.SignalProcessor(
        operations=[sp.IFFT(),sp.apodization.Gaussian(FWHM="100 Hz"),sp.FFT()]
    )

.. skip: next

.. plot::
    :context: close-figs

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.5), subplot_kw={"projection": "csdm"})
    ax[0].plot(processor.apply_operations(dataset=sim.methods[0].simulation))
    ax[0].set_title("Single-Quantum Spectrum All Transitions")
    ax[0].grid()
    ax[0].invert_xaxis()  # reverse x-axis
    ax[1].plot(processor.apply_operations(dataset=sim.methods[1].simulation))
    ax[1].plot(processor.apply_operations(dataset=sim.methods[2].simulation))
    ax[1].set_title("Single-Quantum Spectrum Single Transitions")
    ax[1].grid()
    ax[1].invert_xaxis()  # reverse x-axis
    plt.tight_layout()
    plt.show()


.. note::

    Whenever the ``D`` attribute is omitted, the SymmetryQuery allows
    transitions with all values of :math:`\text{d}_I`. On the other hand, 
    whenever the ``P`` attribute is omitted it takes on a default value 
    of :math:`\text{p}_I = 0`, except when it is omitted from a ``ch1`` SymmetryQuery, in which case it defaults to :math:`\text{p}_I = -1`.

In the example above, we create the SymmetryQuery and TransitionQuery objects
using Python dictionaries.  To do this, the dictionary must use the object's attribute names as the key strings.  


In a notable case, that :math:`\text{d}_I = 0` for all symmetric 
:math:`(m \rightarrow - m)` transitions, is particularly useful for 
half-integer quadrupolar nuclei, as these transitions are unaffected by 
the first-order quadrupolar coupling frequency contribution.  Thus, 
:math:`\ketbra{-\tfrac{1}{2}}{\tfrac{1}{2}}`, the so-called "central transition"
of a quadrupolar nucleus, and the symmetric triple quantum transition 
:math:`\ketbra{-\tfrac{3}{2}}{\tfrac{3}{2}}` are selected in the two-dimensional
custom Method shown below.

.. skip: next

.. plot::
    :context: close-figs

    mqmas = Method(
        channels=["87Rb"],
        magnetic_flux_density=9.4,
        rotor_frequency=10000,
        spectral_dimensions=[
            SpectralDimension(
                count=128,
                spectral_width=6e3,  # in Hz
                reference_offset=-9e3,  # in Hz
                label="3Q resonances",
                events=[
                    SpectralEvent(transition_query=[{"ch1": {"P": [-3], "D": [0]}}])
                ]
            ),
            SpectralDimension(
                count=256,
                spectral_width=6e3,  # in Hz
                reference_offset=-5e3,  # in Hz
                label="1Q resonances",
                events=[
                    SpectralEvent(transition_query=[{"ch1": {"P":[-1], "D": [0]}}])
                ]
            )
        ],
    )

    site1 = Site(
      isotope="87Rb",
      isotropic_chemical_shift=-27.4,  # ppm
      quadrupolar={"Cq":1.68e6, "eta":0.2}  # Cq in Hz
      )
    site2 = Site(
      isotope="87Rb",
      isotropic_chemical_shift=-28.5,  # ppm
      quadrupolar={"Cq":1.94e6, "eta":1}  # Cq in Hz
      )
    site3 = Site(
      isotope="87Rb",
      isotropic_chemical_shift=-31.3,  # ppm
      quadrupolar={"Cq":1.72e6, "eta":0.5}  # Cq in Hz
      )

    # No Couplings, so create a separate SpinSystem for each site.
    sites = [site1, site2, site3]
    spin_systems = [SpinSystem(sites=[s]) for s in sites]

    sim = Simulator(spin_systems=spin_systems, methods=[mqmas])
    sim.run()

    # Apply Gaussian line broadening along both dimensions
    processor = sp.SignalProcessor(
        operations=[
            sp.IFFT(dim_index=(0, 1)),
            sp.apodization.Gaussian(FWHM="0.08 kHz", dim_index=0),
            sp.apodization.Gaussian(FWHM="0.22 kHz", dim_index=1),
            sp.FFT(dim_index=(0, 1)),
        ]
    )
    data = processor.apply_operations(dataset=sim.methods[0].simulation)

    plt.figure(figsize=(6, 4))
    ax = plt.subplot(projection="csdm")
    cb = ax.imshow(data.real / data.real.max(), aspect="auto", cmap="gist_ncar_r")
    plt.colorbar(cb)
    ax.invert_xaxis()
    ax.invert_yaxis()
    plt.tight_layout()
    plt.show()


.. warning::
  This custom Method, as well as the built-in Multi-Quantum VAS methods, 
  assumes uniform excitation and mixing of the multiple-quantum transition.
  In an experimental MQ-MAS measurement both excitation and mixing efficiencies
  are dependent on the ratio of the quadrupolar coupling constant to the
  rf field strength.  Therefore, the relative integrated intensities of this
  simulation may not agree with experiment.


For 3Q-MAS experiments, the spectrum is often sheared and scaled 
to make the vertical dimension the purely isotropic dimension. This can 
be accomplished with an affine matrix added to the method. Let’s re-make 
our 3Q-MAS method with this affine matrix.

.. skip: next

.. plot::
    :context: close-figs

    sheared_mqmas = Method(
        channels=["87Rb"],
        magnetic_flux_density=9.4,
        rotor_frequency=10000,
        spectral_dimensions=[
            SpectralDimension(
                count=128,
                spectral_width=6e3,  # in Hz
                reference_offset=-9e3,  # in Hz
                label="Isotropic dimension",
                events=[
                    SpectralEvent(transition_query=[{"ch1": {"P": [-3], "D": [0]}}])
                ]
            ),
            SpectralDimension(
                count=256,
                spectral_width=6e3,  # in Hz
                reference_offset=-5e3,  # in Hz
                label="MAS dimension",
                events=[
                    SpectralEvent(transition_query=[{"ch1": {"P":[-1], "D": [0]}}])
                ]
            )
        ],
        affine_matrix=[[9/16, 7/16], [0, 1]]
    )

    sim = Simulator(spin_systems=spin_systems, methods=[sheared_mqmas])
    sim.run()

    # Apply Gaussian line broadening along both dimensions
    processor = sp.SignalProcessor(
        operations=[
            sp.IFFT(dim_index=(0, 1)),
            sp.apodization.Gaussian(FWHM="0.08 kHz", dim_index=0),
            sp.apodization.Gaussian(FWHM="0.22 kHz", dim_index=1),
            sp.FFT(dim_index=(0, 1)),
        ]
    )
    data = processor.apply_operations(dataset=sim.methods[0].simulation)

    plt.figure(figsize=(6, 4))
    ax = plt.subplot(projection="csdm")
    cb = ax.imshow(data.real / data.real.max(), aspect="auto", cmap="gist_ncar_r")
    plt.colorbar(cb)
    ax.invert_xaxis()
    ax.invert_yaxis()
    plt.tight_layout()
    plt.show()

See the `"Symmetry Pathways in Solid-State NMR" paper
<https://doi.org/10.1016/j.pnmrs.2010.11.003>`_  for more details on affine transforations.

You may have noticed that the ``transition_queries`` attribute of SpectralEvent 
holds a list of TransitionQuery objects.   Each TransitionQuery in the list 
applies to the full set of transitions present in the spin system.  It is the 
union of these transition subsets that becomes the final set of selected 
transitions during the SpectralEvent. Consider the case of the two satellite 
transitions closest to one side of the central transition of a quadrupolar 
nucleus.  The :math:`\text{p}_I` and :math:`\text{d}_I` values for these 
two transitions are

- :math:`|-1/2\rangle\rightarrow|-3/2\rangle \,\,\,\,\,\text{is}\,\,\,\left(\text{p}_I, \text{d}_I\right)=(-1,2),`
- :math:`|-3/2\rangle\rightarrow|-5/2\rangle \,\,\,\,\,\text{is}\,\,\,\left(\text{p}_I, \text{d}_I\right)=(-1,4).`

These two transitions can be selected using the code below.

.. plot::
    :context: close-figs

    from mrsimulator.method import SpectralEvent

    event = SpectralEvent(
        transition_query=[
            {"ch1": {"P": [-1], "D": [2]}},
            {"ch1": {"P": [-1], "D": [4]}},
        ]
    )



We have seen how a Method object can select between different coherences by using
SpectralDimension and SpectralEvents. By adding a MixingEvent, we can selectively simulate
frequencies from specific transition pathways. Below we construct a deuterium spin system
and two Method objects to simulate a Hahn and Solid Echo experiment.

Hahn Echo
"""""""""

The Hahn Echo experiment observes the transition frequencies from the following
:math:`\mathbb{p}` transition symmetry pathways (a.k.a coherence transfer pathways).

.. math::

    \mathbb{p}: 0 \xrightarrow[]{\frac{\pi}{2}} +1 \xrightarrow[]{\pi} -1

This pathway selectively refocuses the :math:`\mathbb{p}` frequency contributions into
an echo while leaving the :math:`\mathbb{d}` contributions free to evolve unaffected by the
:math:`\pi` pulse.
Below is a diagram representing the different energy level transitions and corresponding
pathways observed by the Hahn Echo experiment.

.. figure:: ../../_static/deuteriumHahnEcho.*
    :alt: Transition symmetry pathways for the Hahn Echo experiment
    :align: center
    :width: 50%

    Energy level transitions and symmetry pathways for the Hahn Echo experiment.

Although a normal experiment would start with a :math:`\frac{\pi}{2}` rotation to transfer the
equilibrium magnetization to a desired symmetry, mrsimulator eliminates the need for this first
rotation by defining the first symmetry as :math:`\mathbb{p} = +1`. Our transition symmetry
pathway now becomes

.. math::

    \mathbb{p}: +1 \xrightarrow[]{\pi} -1

Below is a method object which simulated the Hahn Echo experiment. The MixingEvent defines the
:math:`\pi` rotation between the two SpectralEvents.

.. plot::
    :context: close-figs

    from mrsimulator.method import MixingEvent
    from pprint import pprint

    hahn_echo = Method(
        channels=["2H"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            SpectralDimension(
                count=512,
                spectral_width=2e4,  # in Hz
                events=[
                    SpectralEvent(fraction=0.5, transition_query=[
                        {"ch1": {"P": [1], "D": [1]}},
                        {"ch1": {"P": [1], "D": [-1]}},
                    ]),
                    MixingEvent(query={"ch1": {"angle": 3.141592, "phase": 0}}),
                    SpectralEvent(fraction=0.5, transition_query=[
                        {"ch1": {"P": [-1], "D": [1]}},
                        {"ch1": {"P": [-1], "D": [-1]}},
                    ])
                ]
            )
        ]
    )

    spin_system = SpinSystem(sites=[deuterium])
    pprint(hahn_echo.get_transition_pathways(spin_system))

.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    [|1.0⟩⟨0.0| ⟶ |-1.0⟩⟨0.0|, weight=(1+0j)
     |0.0⟩⟨-1.0| ⟶ |0.0⟩⟨1.0|, weight=(1+0j)]



Solid Echo
""""""""""

The Solid Echo experiment selectively refocuses the the :math:`\mathbb{d}` frequency contributions
into an echo using a :math:`\frac{\pi}{2}` rotation while keeping the :math:`\mathbb{p}` pathway
constant.
Below is a diagram representing the different energy level transitions and corresponding
pathways observed by the Solid Echo experiment.

.. figure:: ../../_static/deuteriumSolidEcho.*
    :alt: Transition symmetry pathways for the Hahn Echo experiment
    :align: center
    :width: 50%

    Energy level transitions and symmetry pathways for the Solid Echo experiment.

.. math::

    \mathbb{p}: -1 \xrightarrow[]{\frac{\pi}{2}} -1

    \mathbb{d}: \pm 1 \xrightarrow[]{\frac{\pi}{2}} \mp 1

Below we construct the Solid Echo method and print out the transition pathways for the
deuterium spin system.

.. plot::
    :context: close-figs

    solid_echo = Method(
        channels=["2H"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            SpectralDimension(
                count=512,
                spectral_width=2e4,  # in Hz
                events=[
                    SpectralEvent(fraction=0.5, transition_query=[
                        {"ch1": {"P": [-1], "D": [1]}},
                        {"ch1": {"P": [-1], "D": [-1]}},
                    ]),
                    MixingEvent(query={"ch1": {"angle": 3.141592 / 2, "phase": 0}}),
                    SpectralEvent(fraction=0.5, transition_query=[
                        {"ch1": {"P": [-1], "D": [1]}},
                        {"ch1": {"P": [-1], "D": [-1]}},
                    ]),
                ]
            )
        ]
    )

    pprint(solid_echo.get_transition_pathways(spin_system))

.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    [|-1.0⟩⟨0.0| ⟶ |0.0⟩⟨1.0|, weight=(0.5+0j)
     |0.0⟩⟨1.0| ⟶ |-1.0⟩⟨0.0|, weight=(0.5+0j)]

.. note::

    Although we explicitly defined the :math:`D` values for each transition query in the
    above method, mrsimulator will expand an undefined :math:`D` to all allowed values.
    The transition queries in the Solid Echo method could have just as easily been defined
    as ``{"ch1": {"P": [-1]}}``.

Now we setup and run the simulation then process and plot the data

.. skip: next

.. plot::
    :context: close-figs
    :caption: Simulated Hahn Echo spectrum (left) and Solid Echo spectrum (right) for the same :math:`2^\text{H}` spin system.

    sim = Simulator()
    sim.spin_systems = [spin_system]
    sim.methods = [hahn_echo, solid_echo]
    sim.run()

    processor = sp.SignalProcessor(
        operations=[
            sp.IFFT(),
            sp.apodization.Gaussian(FWHM="100 Hz"),
            sp.FFT(),
        ]
    )
    hahn_data = processor.apply_operations(dataset=sim.methods[0].simulation)
    solid_data = processor.apply_operations(dataset=sim.methods[1].simulation)

    fig, ax = plt.subplots(1, 2, subplot_kw={"projection": "csdm"}, figsize=[8.5, 3])
    ax[0].plot(hahn_data.real, color="black", linewidth=1)
    ax[0].invert_xaxis()
    ax[1].plot(solid_data.real, color="black", linewidth=1)
    ax[1].invert_xaxis()
    plt.tight_layout()
    plt.show()

    
Multi-Spin Transition Queries
-----------------------------

When there is more than one site in a spin system, things get a little more 
complicated. 

Single-Spin Single-Quantum Transitions
''''''''''''''''''''''''''''''''''''''

Consider the case of three weakly coupled proton sites.  Here, the
selection rule for observable transitions is

.. math::
    \left.
    \begin{array}{ll}
    \text{p}_A = - 1 \mbox{  while  }  \text{p}_M = 0, \text{p}_X = 0 \\
    \text{p}_M = - 1 \mbox{  while  }  \text{p}_A = 0, \text{p}_X = 0 \\
    \text{p}_X = - 1 \mbox{  while  }  \text{p}_A = 0, \text{p}_M = 0 \\
    \end{array}
    \right\}
    \text{ Detection Selection Rules.}

These corresponds to the *single-spin 
single-quantum transitions* labeled :math:`\hat{A}_1`, 
:math:`\hat{A}_2`, :math:`\hat{A}_3`, :math:`\hat{A}_4`, :math:`\hat{M}_1`, 
:math:`\hat{M}_2`, :math:`\hat{M}_3`,:math:`\hat{M}_4`, :math:`\hat{X}_1`,
math:`\hat{X}_2`, math:`\hat{X}_3`, and :math:`\hat{X}_4` 
in the energy level diagram below.

.. figure:: ../../_static/ThreeCoupledSpinsEnergy.*
    :width: 600
    :alt: figure
    :align: center

Keep in mind that the Method object does not know, in advance, the 
number of sites in a spin system.   

The TransitionQuery for selecting these 12 *single-spin single-quantum* transitions
is given in the code below.

.. skip: next

.. plot::
    :context: close-figs
    
    site_A = Site(isotope="1H", isotropic_chemical_shift=0.5)
    site_M = Site(isotope="1H", isotropic_chemical_shift=2.5)
    site_X = Site(isotope="1H", isotropic_chemical_shift=4.5)
    sites = [site_A,site_M,site_X]
    coupling_AM = Coupling(site_index=[0, 1], isotropic_j=12)
    coupling_AX = Coupling(site_index=[0, 2], isotropic_j=12)
    coupling_MX = Coupling(site_index=[1, 2], isotropic_j=12)
    couplings = [coupling_AM, coupling_AX, coupling_MX]
    system = SpinSystem(sites=sites, couplings=couplings)

    method = Method(
        channels=["1H"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            SpectralDimension(
                count=16000,
                spectral_width=1800,  # in Hz
                reference_offset=1000,  # in Hz
                label="$^{1}$H frequency",
                events=[
                  {
                  "fraction":1,
                  "transition_query":[{"ch1":{"P":[-1]}}]
                  }
                ]
            )
        ]
    )

    sim = Simulator(spin_systems = [system],methods=[method])
    sim.run()

    processor = sp.SignalProcessor(
        operations=[sp.IFFT(),sp.apodization.Exponential(FWHM="1 Hz"),sp.FFT()]
    )

    plt.figure(figsize=(10, 3))  # set the figure size
    ax = plt.subplot(projection="csdm")
    ax.plot(processor.apply_operations(dataset=sim.methods[0].simulation))
    ax.invert_xaxis()  # reverse x-axis
    plt.tight_layout()
    plt.grid()
    plt.show()

The assignment of transitions in the spectrum above are, from left to right, are 
:math:`\hat{X}_4, (\hat{X}_3, \hat{X}_2)`, and :math:`\hat{X}_1` centered at 
4.5 ppm, :math:`\hat{M}_4, (\hat{M}_3, \hat{M}_2)`, and :math:`\hat{M}_1` 
centered at 2.5 ppm, and :math:`\hat{A}_4, (\hat{A}_3, \hat{A}_2)`, and 
:math:`\hat{A}_1` centered at 0.5 ppm.

To understand how the TransitionQuery works in this case, it is important to 
realize that all Sites having same isotope are  "indistinguishable" to a 
TransitionQuery object.  Recall that ``ch1`` is associated with the first 
isotope in the list of isotope strings assigned to the Method attribute 
``channels``.   When the TransitionQuery above is combined with the SpinSystem 
object with three :math:`^1\text{H}` Sites, it must first expand its 
SymmetryQuery into an intermediate set of spin-system-specifc 
symmetry queries, illustrated by each row in the table below.

.. list-table:: 
   :widths: 25 25 25 25
   :header-rows: 1

   * - Transitions
     - :math:`\text{p}_A`
     - :math:`\text{p}_M`
     - :math:`\text{p}_X`
   * - :math:`\hat{A}_1, \hat{A}_2, \hat{A}_3, \hat{A}_4`
     - –1
     - 0
     - 0
   * - :math:`\hat{M}_1, \hat{M}_2, \hat{M}_3, \hat{M}_4`
     - 0
     - –1
     - 0
   * - :math:`\hat{X}_1, \hat{X}_2, \hat{X}_3, \hat{X}_4`
     - 0
     - 0
     - –1

The intermediate spin-system-specifc symmetry query in each row is used to 
select a subset of transitions from the full set of transitions.  The
final set of selected transitions is obtained from the union of transition 
subsets from each spin-system-specifc symmetry query.

The :py:meth:`~mrsimulator.Method.get_transition_pathways` function will allow
you to inspect the transitions selected by the TransitionQuery objects in the
SpectralEvent in terms of the initial and final Zeeman eigenstate quantum
numbers.

.. skip: next

.. plot::
    :context: close-figs

    from pprint import pprint
    pprint(method.get_transition_pathways(system))


.. figure:: ../../_static/method_user_doc_pprint_output1.*
    :width: 450
    :alt: figure
    :align: center


To further illustrate how the TransitionQuery and SymmetryQuery objects 
works in a multi-site spin system, let's examine a few more examples in 
the case of three weakly coupled proton sites.   

Two-Spin Double-Quantum Transitions
'''''''''''''''''''''''''''''''''''

In this spin system there are six *two-spin double-quantum transitions* where 
:math:`\text{p}_{AMX} = \text{p}_{A} + \text{p}_{M} + \text{p}_{X} = -2` and
another six *two-spin double-quantum transitions* where 
:math:`\text{p}_{AMX} = \text{p}_{A} + \text{p}_{M} + \text{p}_{X} = +2`.  The
:math:`\text{p}_{AMX} = -2` transitions are illustrated in the energy-level diagram 
below.

.. figure:: ../../_static/ThreeCoupledSpinsDoubleQuantum.*
    :width: 600
    :alt: figure
    :align: center


The code below will select the six *two-spin double-quantum transitions* where 
:math:`\text{p}_{AMX} = -2`.

.. skip: next

.. plot::
    :context: close-figs

    method = Method(
        channels=["1H"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            SpectralDimension(
                count=16000,
                spectral_width=2000,  # in Hz
                reference_offset=2000,  # in Hz
                label="$^{1}$H frequency",
                events=[
                  {
                  "fraction":1,
                  "transition_query":[{"ch1":{"P":[-1,-1]}}]
                  }
                ]
            )
        ]
    )

    sim = Simulator(spin_systems = [system],methods=[method])
    sim.run()

    plt.figure(figsize=(10, 3))  # set the figure size
    ax = plt.subplot(projection="csdm")
    ax.plot(processor.apply_operations(dataset=sim.methods[0].simulation))
    ax.invert_xaxis()  # reverse x-axis
    plt.tight_layout()
    plt.grid()
    plt.show()

The assignment of transitions in the spectrum above are, from left to right, are 
:math:`\hat{D}_{2,MX}`, :math:`\hat{D}_{1,MX}`, :math:`\hat{D}_{2,AX}`, 
:math:`\hat{D}_{1,AX}`, :math:`\hat{D}_{2,AM}`, and
:math:`\hat{D}_{1,AM}`, 


As before, when this generic TransitionQuery is combined with the three-site 
SpinSystem object, the SymmetryQuery is expanded into an intermediate set of 
spin-system-specifc symmetry queries illustrated in the table below.

.. list-table:: 
   :widths: 25 25 25 25
   :header-rows: 1

   * - Transitions
     - :math:`\text{p}_A`
     - :math:`\text{p}_M`
     - :math:`\text{p}_X`
   * - :math:`\hat{D}_{1,AM}, \hat{D}_{2,AM}`
     - –1
     - –1
     - 0
   * - :math:`\hat{D}_{1,MX}, \hat{D}_{2,MX}`
     - 0
     - –1
     - –1
   * - :math:`\hat{D}_{1,AX}, \hat{D}_{2,AX}`
     - –1
     - 0
     - –1

Again, the intermediate spin-system-specifc symmetry query in each row is used to 
select a subset of transitions from the full set of transitions.  The
final set of selected transitions is obtained from the union of transition 
subsets from each spin-system-specifc symmetry query.

.. skip: next

.. plot::
    :context: close-figs

    from pprint import pprint
    pprint(method.get_transition_pathways(system))


.. figure:: ../../_static/method_user_doc_pprint_output2.*
    :width: 450
    :alt: figure
    :align: center



Three-Spin Single-Quantum Transitions
'''''''''''''''''''''''''''''''''''''

Another interesting example in this spin system with three weakly coupled 
proton sites are the three *three-spin single-quantum transitions* having 
:math:`\text{p}_{AMX} = \text{p}_{A} + \text{p}_{M} + \text{p}_{X} = -1` and the
three *three-spin single-quantum transitions* having 
:math:`\text{p}_{AMX} = \text{p}_{A} + \text{p}_{M} + \text{p}_{X} = +1`

The three *three-spin single-quantum transitions* having 
:math:`\text{p}_{AMX} = -1` are illustrated in the energy level diagram below.

.. figure:: ../../_static/ThreeCoupledSpinsSingleQuantum.*
    :width: 600
    :alt: figure
    :align: center

The code below will select these *three-spin single-quantum transitions*.

.. skip: next

.. plot::
    :context: close-figs

    method = Method(
        channels=["1H"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            SpectralDimension(
                count=16000,
                spectral_width=4000,  # in Hz
                reference_offset=1000,  # in Hz
                label="$^{1}$H frequency",
                events=[
                  {
                  "fraction":1,
                  "transition_query":[{"ch1":{"P":[-1,-1,+1]}}]
                  }
                ]
            )
        ]
    )

    sim = Simulator(spin_systems = [system],methods=[method])
    sim.run()

    plt.figure(figsize=(10, 3))  # set the figure size
    ax = plt.subplot(projection="csdm")
    ax.plot(processor.apply_operations(dataset=sim.methods[0].simulation))
    ax.invert_xaxis()  # reverse x-axis
    plt.tight_layout()
    plt.grid()
    plt.show()

The assignment of transitions in the spectrum above are, from left to right, are 
:math:`\hat{S}_{3,AMX}`, :math:`\hat{S}_{2,AMX}`, and :math:`\hat{S}_{1,AMX}`

Again, combined with the three-site SpinSystem object, the SymmetryQuery is 
expanded into the set of spin-system-specifc symmetry queries illustrated 
in the table below.

.. list-table:: 
   :widths: 25 25 25 25
   :header-rows: 1

   * - Transitions
     - :math:`\text{p}_A`
     - :math:`\text{p}_M`
     - :math:`\text{p}_X`
   * - :math:`\hat{S}_{1,AMX}`
     - –1
     - +1
     - –1
   * - :math:`\hat{S}_{2,AMX}`
     - –1
     - –1
     - +1
   * - :math:`\hat{S}_{3,AMX}`
     - +1
     - –1
     - –1



.. skip: next

.. plot::
    :context: close-figs

    from pprint import pprint
    pprint(method.get_transition_pathways(system))


.. figure:: ../../_static/method_user_doc_pprint_output3.*
    :width: 450
    :alt: figure
    :align: center


As you can surmise from the examples, the attributes of 
SymmetryQuery, ``P`` and ``D``, hold a list of single-spin transition 
symmetry function values, and the length of the list is the desired number 
of spins that are involved in the transition.

Heteronuclear multiple-spin transitions
'''''''''''''''''''''''''''''''''''''''

How does ``D`` fit into the multi-site SymmetryQuery story? Consider the 
case of two coupled hydrogen, except we replace one of the :math:`^1H` with
:math:`^2H`.  Let's focus on the single-spin single-quantum transitions, shown below as :math:`\hat{A}_{1\pm}` and :math:`\hat{A}_{2\pm}` on the left, and the two-spin triple-quantum transition, shown below as  :math:`\hat{T}_{AX}` on
the right.

.. figure:: ../../_static/Spin1SpinHalfCouple.*
    :width: 900
    :alt: figure
    :align: center


.. skip: next

.. plot::
    :context: close-figs

    import numpy as np

    site_A = Site(isotope="2H", isotropic_chemical_shift=0.5, 
      quadrupolar={
          "Cq":100000,  # in Hz
          "eta":0.2,
          "alpha":5 * np.pi / 180,
          "beta":np.pi / 2,
          "gamma":70 * np.pi / 180}
          )
    site_X = Site(isotope="1H", isotropic_chemical_shift=4.5)
    sites = [site_A,site_X]
    coupling_AX = Coupling(site_index=[0, 1], dipolar={"D":-20000})
    couplings = [coupling_AX]
    system = SpinSystem(sites=sites, couplings=couplings)

    methodAll1Q = Method(
      channels=["2H","1H"],
      magnetic_flux_density=9.4,  # in T
      spectral_dimensions=[
          SpectralDimension(
              count=16000,
              spectral_width=200000,  # in Hz
              reference_offset=0,  # in Hz
              label="$^{2}$H frequency",
              events=[
                {
                "fraction":1,
                "transition_query":[{"ch1":{"P":[-1]}}]
                }
              ]
          )
      ]
    )

    methodHalf1Q = Method(
      channels=["2H","1H"],
      magnetic_flux_density=9.4,  # in T
      spectral_dimensions=[
          SpectralDimension(
              count=16000,
              spectral_width=200000,  # in Hz
              reference_offset=0,  # in Hz
              label="$^{2}$H frequency",
              events=[
                {
                "fraction":1,
                "transition_query":[{"ch1":{"P":[-1],"D":[-1]}}]
                }
              ]
          )
      ]
    )

    method3Q = Method(
        channels=["2H","1H"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            SpectralDimension(
                count=16000,
                spectral_width=10000,  # in Hz
                reference_offset=5000,  # in Hz
                label="$^{2}$H frequency",
                events=[
                    {
                        "fraction":1,
                        "transition_query":[{
                            "ch1":{"P":[-2]},"ch2":{"P":[-1]}
                        }]
                    }
                ]
            )
        ]
    )
    processor = sp.SignalProcessor(
      operations=[sp.IFFT(),sp.apodization.Gaussian(FWHM="100 Hz"),sp.FFT()]
    )

    sim = Simulator(spin_systems=[system],methods=[methodAll1Q,methodHalf1Q,method3Q])
    sim.config.integration_volume = "hemisphere"
    sim.run()

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.5), subplot_kw={"projection": "csdm"})
    ax[0].plot(processor.apply_operations(dataset=sim.methods[0].simulation))
    ax[0].set_title("Full Single-Quantum Spectrum")
    ax[0].grid()
    ax[0].invert_xaxis()  # reverse x-axis
    ax[1].plot(processor.apply_operations(dataset=sim.methods[1].simulation))
    ax[1].set_title("Half Single-Quantum Spectrum")
    ax[1].grid()
    ax[1].invert_xaxis()  # reverse x-axis
    plt.tight_layout()
    plt.show()


The deuterium spectrum of a static-polycrystalline sample is shown on the left is for all single-spin single-quantum transitions on deuterium, :math:`\hat{A}_{1\pm}` and :math:`\hat{A}_{2\pm}`.  The spectrum on the right is for half of the single-spin single-quantum transitions on deuterium: :math:`\hat{A}_{1-}` and :math:`\hat{A}_{2-}`.

.. skip: next

.. plot::
    :context: close-figs

    plt.figure(figsize=(10, 3))  # set the figure size
    ax = plt.subplot(projection="csdm")
    ax.set_title("Heteronuclear Two-Spin ($^2$H-$^1$H) Triple-Quantum Spectrum")
    ax.plot(processor.apply_operations(dataset=sim.methods[2].simulation))
    plt.tight_layout()
    plt.grid()
    plt.show()

.. list-table:: 
   :widths: 25 25 25 25
   :header-rows: 1

   * - Transitions
     - :math:`\text{p}_A`
     - :math:`\text{d}_A`
     - :math:`\text{p}_X`
   * - :math:`\hat{T}_{AX}`
     - –2
     - 0
     - –1

The single transition in the heteronuclear two-spin (:math:`^2\text{H}`-:math:`^1\text{H}`) triple-quantum spectrum is unaffected by the dipolar and quadrupolar frequency anisotropies.

Mixing Events
-------------

Default Mixing between Events
'''''''''''''''''''''''''''''

.. plot::
    :context: close-figs


    from mrsimulator.method import MixingEvent

    MixingEvent(query="TotalMixing"),


No Mixing Event
'''''''''''''''

.. plot::
    :context: close-figs


    MixingEvent(query="NoMixing"),

Rotational Query
''''''''''''''''

Mixing events are used to transfer (permute) among transitions and populations,
e.g., :math:`\pi/2` or :math:`\pi` rotations between consecutive spectral or
delay events.  For a rotation in a mixing event, the efficiency
associated with the coherence transfer from 

.. math::
    :label: transition

    \ketbra{I, m_f}{I, m_i} \stackrel{\theta_\phi}{\longrightarrow} a(I,\theta,\phi) \ketbra{I,m_f'}{I,m_i'}

is 

.. math::
    :label: rotation

     a(I,\theta,\phi) = d_{m_f',m_f}^{(I)}(\theta)d_{m_i',m_i}^{(I)}(\theta)e^
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

.. plot::
    :context: close-figs


    import numpy as np
    from mrsimulator.method.query import RotationalQuery
    rot_query = RotationalQuery(angle=np.pi/2, phase=0)
    rot_mixing = MixingEvent(query={"ch1": rot_query})

It is through :py:meth:`~mrsimulator.method.query.MixingQuery` and 
:py:meth:`~mrsimulator.method.query.TransitionQuery` 
objects that the desired transition pathways are selected and undesired transition 
pathways are eliminated.


Frequency Contributions
-----------------------


Transition and Symmetry Pathways
--------------------------------

The number of possible transition pathways for a spin system depends on the
number of energy eigenstates and the number of spectral and delay
events in a method. 



SpectralDimension
-----------------

Mrsimulator allows users to create custom methods and simulate the NMR spectrum.
At the top level, a :ref:`method_api` object is no different than the pre-built
methods provided within the ``mrsimulator.method.lib`` module.

A generic setup for a custom method (similar to the stock method) follows,

.. plot::
    :context: close-figs


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

.. plot::
    :context: close-figs


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

.. plot::
    :context: close-figs


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

.. plot::
    :context: close-figs

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
attributes that have already been discussed in earlier sections of the documentation. 
Notably, ``origin_offset`` and ``reference_offset`` are important for converting
the frequency coordinate into a dimensionless frequency ratio coordinate. For
spectra where all the spectral dimensions are associated with single-quantum
transitions on a single isotope, the convention for defining ``origin_offset`` 
and ``reference_offset`` is well established;
the ``origin_offset``, :math:`o_k`, is interpreted as the NMR spectrometer
frequency and  the ``reference_offset``, :math:`b_k`, as the reference
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

.. plot::
    :context: close-figs

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

.. plot::
    :context: close-figs

    SpectralDimension(
        events=[
            SpectralEvent(fraction=0.5, transition_query=[{"ch1": {"P": [1]}}]),
            MixingEvent(query={"ch1": {"angle": 3.14159, "phase": 0}}),
            SpectralEvent(fraction=0.5, transition_query=[{"ch1": {"P": [-1]}}]),
        ]
    )

**A one-dimensional solid echo**

:math:`\mathbb{p}: -1 \xrightarrow[]{\frac{\pi}{2}} -1`

.. plot::
    :context: close-figs

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
.. list-table:: The attributes of a Method object
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
