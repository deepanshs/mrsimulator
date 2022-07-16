
.. _method_documentation:

======
Method
======

While **mrsimulator**'s organization of the :ref:`spin_sys_api` object and its
composite objects, :ref:`site_api` and :ref:`coupling_api`, are easily
understood by anyone familiar with the underlying physical concepts, the
organization of the :ref:`method_api` object in **mrsimulator** and its related
composite objects require a more detailed explanation of their design. This
section assumes that you are already familiar with the topics covered in the
Introduction sections :ref:`getting_started`,
:ref:`introduction_isotopomers_example`, and :ref:`fitting_example`.

.. note::

    Before writing your own custom Method, check if any of our pre-built methods in the :ref:`methods_library_documentation` can serve your needs.


Overview
--------

An experimental NMR method involves a sequence of rf pulses, free evolution
periods, and sample motion. The Method object in **mrsimulator** models the
spectrum from an NMR pulse sequence. The Method object is designed to be
versatile in its ability to model spectra from various multi-pulse NMR methods
using concepts from the `symmetry pathway approach
<https://doi.org/10.1016/j.pnmrs.2010.11.003>`_ where a pulse sequence is
understood in terms of a set of desired (and undesired)  *transition pathways*.
Each transition pathway is associated with a single resonance in a
multi-dimensional NMR spectrum. The transition pathway signal encodes
information about the spin system interactions in its amplitude and correlated
frequencies. Consider the illustration of a 2D pulse sequence shown below, where
a desired signal for the method is associated with a particular transition
pathway, :math:`{\hat{A} \rightarrow \hat{B} \rightarrow \hat{C} \rightarrow
\hat{D} \rightarrow \hat {E} \rightarrow \hat{F}}`.

.. figure:: ../../_static/TransitionPathway.*
    :width: 600
    :alt: figure
    :align: center

    A illustration of an two-dimensional NMR pulse sequence leading up to the
    acqusition of the signal from a transition pathway.


Here, the first spectral dimension, i.e., the Fourier transform of the
transition pathway signal as a function of :math:`t_1`, derives its *average
frequency*, :math:`\overline{\Omega}_1`, from a weighted average of the
:math:`\hat{A}`, :math:`\hat{B}`, and :math:`\hat{C}` transition frequencies.
The second spectral dimension, i.e., the FT with respect to :math:`t_2`, derives
its average frequency, :math:`\overline{\Omega}_2`, from a weighted average of
the :math:`\hat{E}`, and :math:`\hat{F}` transition frequencies. Much of the
experimental design and implementation of an NMR method is in identifying the
desired transition pathways and finding ways to acquire their signals while
eliminating all undesired transition pathway signals.

While NMR measurements occur in the time domain, **mrsimulator** simulates the
corresponding multi-dimensional spectra directly in the frequency domain. The
Method object in **mrsimulator** needs only a few details of the NMR pulse
sequence to generate the spectrum. It mimics the result of the pulse sequence
given the desired transition pathways and their complex amplitudes and average
frequencies in each spectroscopic dimension of the dataset. To this end, a
Method object is organized according to the UML diagram below.

.. figure:: ../../_static/MethodUML.*
    :width: 700
    :alt: figure
    :align: center

    Unified Modeling Language class diagram of the Method object in mrsimulator.

.. note::

 In UML (Unified Modeling Language) diagrams, each class is represented with a
 box that contains two compartments. The top compartment has the class's name,
 and the bottom compartment contains the class's attributes. Default attribute
 values are shown as assignments. A composition is depicted as a binary
 association decorated with a filled black diamond. Inheritance is shown as a
 line with a hollow triangle as an arrowhead.

At the heart of a Method object, assigned to its attribute
``spectral_dimensions``, is an ordered list of :ref:`spectral_dim_api` objects
in the same order as the time evolution dimensions of the experimental NMR
sequence. In each SpectralDimension object, an ordered list of :ref:`event_api`
objects assigned to the attribute ``events``; Event objects are divided
into three types: (1) :py:meth:`~mrsimulator.method.SpectralEvent`, (2)
:py:meth:`~mrsimulator.method.DelayEvent`, and (3)
:py:meth:`~mrsimulator.method.MixingEvent`.  This ordered list of Event objects
is used to select the desired transition pathways and determine their average
frequency and complex amplitude in the SpectralDimension.

.. warning::

  DelayEvent objects are not available in version 0.7 of **mrsimulator**.

SpectralEvent and DelayEvent objects define which transitions are
observed during the event and under which transition-dependent frequency
contributions they evolve. No coherence transfer among transitions or
populations occurs in a spectral or delay event. The transition-dependent
frequency contributions during an Event are selected from a list of
:ref:`enumeration literals<freq_contrib_api>` and placed in the ``freq_contrib``
attribute of the event. If ``freq_contrib`` is left unspecified, i.e., the
value of ``freq_contrib`` is set to ``None``, a default list holding the
enumeration literals for *all* contributions is generated for the event.

.. note::

  All frequency contributions from direct and indirect spin-spin couplings are
  calculated in the weak-coupling limit in **mrsimulator**.

Additionally, the user can affect transition frequencies during a spectral or
delay event by changing other measurement attributes: ``rotor_frequency``,
``rotor_angle``, and ``magnetic_flux_density``. If left unspecified, these
attributes default to the values of the identically named global attributes in
the Method object. SpectralEvent objects use the ``fraction`` attribute to
calculate the weighted average frequency during the spectral dimension for each
selected transition pathway.

Inside SpectralEvent and DelayEvent objects, is a list of
:py:meth:`~mrsimulator.method.query.TransitionQuery` objects (*vide infra*)
which determine which transitions are observed during the event. Method
objects in **mrsimulator** are general-purpose because they are designed for an
arbitrary spin system, i.e., a method does not know the spin system in advance.
When designing a Method object, you cannot identify and select a transition
through its initial and final eigenstate quantum numbers. Transition selection
is done through TransitionQuery and
:py:meth:`~mrsimulator.method.query.SymmetryQuery` objects during individual
spectral or delay events. TransitionQuery objects can hold a
SymmetryQuery object in the attributes ``ch1``, ``ch2``, or ``ch3``, which act on
specific isotopes defined by the ``channels`` attribute in Method. It is
only during a simulation that the Method object uses its TransitionQuery
objects to determine the selected transition pathways for a given SpinSystem
object by the initial and final eigenstate quantum numbers of each transition.

Between adjacent SpectralEvent or DelayEvent objects, **mrsimulator** defaults
to *total mixing*, i.e., connecting all selected transitions in the two adjacent
spectral or delay events. This default behavior can be overridden by placing an
explicit MixingEvent object between such events. Inside MixingEvent
objects is a :py:meth:`~mrsimulator.method.query.MixingQuery` object, which
determines the coherence transfer amplitude between transitions. A
MixingQuery object holds
:py:meth:`~mrsimulator.method.query.RotationQuery` objects acting on specific
isotopes in the spin system. As before, the isotope upon which the
RotationQuery objects act is determined by the ``channels`` attribute in the
Method object.

In this guide to designing custom Method objects, we begin with a brief review
of the relevant *Symmetry Pathway* concepts employed in **mrsimulator**. This
review is necessary for understanding (1) how transitions are selected during
spectral and delay events and (2) how average signal frequencies and amplitudes
in each spectral dimension are determined. We outline the procedures for
designing and creating TransitionQuery and MixingQuery for single- and
multi-spin transitions and how to use them to select the transition pathways
with the desired frequency and amplitudes in each SpectralDimension of your
custom Method object. In multi-dimensional spectra, we illustrate how the
desired frequency correlation can sometimes be achieved by using an appropriate
affine transformation. We also examine how changing the frequency contributions
in SpectralEvent of DelayEvent objects can be used to obtain the desired
frequency and amplitude behavior. The ability to select :ref:`frequency
contributions<freq_contrib_api>` can often reduce the number of events needed in
the design of your custom Method object.

Theoretical Background
----------------------

Before giving details on how to create a custom Method object, we review a
few key concepts about spin transitions and *transition symmetry functions*.

The number of quantized energy eigenstates for :math:`N` coupled nuclei is

.. math::

    \Upsilon_{\left\{ I_1, I_2, \ldots, I_N \right\}} = \prod_{u=1}^N (2 I_u+1),

where :math:`I_u` is the total spin angular momentum quantum number of the
:math:`u\text{th}` nucleus. The transition from quantized energy level
:math:`E_i` to :math:`E_j` is one of

.. math::

    \mathcal{N}_{\left\{ I_1, I_2, \ldots, I_N \right\}} = \frac{\Upsilon_{\left\{ I_1, I_2, \ldots, I_N \right\}}!}{(\Upsilon_{\left\{ I_1, I_2, \ldots, I_N \right\}}-2)!}

possible transitions between the :math:`\Upsilon` levels.   Here we count
:math:`i  \rightarrow  j` and :math:`j  \rightarrow  i` as different
transitions.  For example, a single spin with angular momentum quantum number :math:`I=3/2`
will have :math:`\Upsilon_{\left\{ 3/2 \right\}} = 2I+1 = 4` energy levels and
:math:`\mathcal{N}_{\left\{ 3/2 \right\}} = 2I(2I+1) = 12` possible NMR
transitions.   A two spin system, with quantum numbers :math:`I = 1/2` and :math:`S = 1/2`,
will have

.. math::
    \Upsilon_{\left\{ 1/2, 1/2 \right\}} = (2I +1) \cdot (2S +1) = 4

energy levels and

.. math::
  \mathcal{N}_{\left\{ 1/2,1/2 \right\}} =
  \frac{[(2I +1) \cdot (2S +1)]!}{((2I +1) \cdot (2S +1)-2)!}
  = \frac{[2 \cdot 2]!}{(2 \cdot 0)!} = 12

possible NMR transitions. We write a transition (coherence) from state :math:`i`
to :math:`j` using the outer product notation :math:`\ketbra{j}{i}`. In
**mrsimulator**, all simulations are performed in the high-field limit and
further, assume that all spin-spin couplings are in the weak limit.

To write a custom Method in **mrsimulator**, you'll need to determine the desired
transition pathways and select the desired transitions during each
SpectralEvent or DelayEvent. Keep in mind, however, that Method
objects are designed without any details of the spin systems upon which they
will act. For example, in the density matrix of a spin system ensemble, one
could easily identify a transition by its row and column indexes. However, those
indexes depend on the spin system and how the spins and their eigenstates have
been assigned to those indexes. Instead, we need a spin-system agnostic approach
for selecting transitions.


.. _transition_query_documentation:

Spin Transition Symmetry Functions
''''''''''''''''''''''''''''''''''

One way you can select a subset of single-spin transitions if you don't
know the energy eigenstate quantum numbers is to request all transitions whose
single-spin transition symmetry function, :math:`\text{p}_I`
is :math:`-1`, i.e.,

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
    in sign, the convention is to only present the :math:`{\text{p}_I = - 1}`
    transition resonances in single-quantum spectra.

By selecting only single-spin transitions with :math:`\text{p}_I = -1`, you get
all the "observed" transitions from the set of all possible transitions.
Similarly, you can use  :math:`\text{p}_I` to select any subset of single-spin
transitions, such as double-quantum :math:`(\text{p}_I = \pm 2)` transitions,
triple-quantum :math:`(\text{p}_I = \pm 3)` transitions, etc.

While specifying :math:`\text{p}_I` alone is not enough to select an individual
single-spin transition, any individual single-spin transition can be
identified by a combination of :math:`\text{p}_I` and the single-spin
transition symmetry function :math:`\text{d}_I`, given by

.. math::

    \text{d}_I(m_i,m_j) =  ~m_j^2 - m_i^2.

You can verify this from the values of :math:`\text
{p}_I` and :math:`\text{d}_I` for all single-spin transitions
for :math:`I=1`, :math:`I=3/2` and :math:`I=5/2` shown below.  Note
that :math:`\text{d}_I = 0` for all transitions in a :math:`I=1/2` nucleus.

.. figure:: ../../_static/SpinOneThreeHalves.*
    :width: 500
    :alt: figure
    :align: center

    Energy level diagrams of a spin :math:`I=1` nucleus  (left) and spin
    :math:`I=3/2` nucleus (right). Arrows beginning at the initial state and end
    at the final state represent transitions.   Transitions are labeled with
    their corresponding :math:`\text{p}_I` and :math:`\text{d}_I` transition
    symmetry function values.

.. figure:: ../../_static/SpinFiveHalf.*
    :width: 550
    :alt: figure
    :align: center

    Energy level diagram of a spin :math:`I=5/2` nucleus. Arrows beginning at the initial state and end at the final state represent transitions.   Transitions are labeled with their corresponding :math:`\text{p}_I` and :math:`\text{d}_I` spin transition symmetry function values.

----

For a summary on spin transition symmetry functions in NMR, click on the disclosure
button below.

.. raw:: html

    <details class="reveal-theory">
    <summary><b>Spin Transition Symmetry Functions</b></summary>

.. include:: spin_trans_symm.rst

.. raw:: html

    </details>

Single-Spin Queries
-------------------

Based on the review above, we now know for the spin :math:`I=1`, the transition
:math:`\ketbra{-1}{0}` can be selected with :math:`{(\text{p}_I,\text{d}_I) =
(-1,1)}`.  In **mrsimulator**, this transition is selected during a
SpectralEvent using the SymmetryQuery and TransitionQuery objects,
as defined in the code below.

.. plot::
    :context: reset

    from mrsimulator.method.query import SymmetryQuery, TransitionQuery
    from mrsimulator.method import SpectralEvent

    symm_query = SymmetryQuery(P=[-1], D=[1])
    trans_query = TransitionQuery(ch1=symm_query)
    spec_event = SpectralEvent(transition_queries=[trans_query])

.. note::
    Python dictionaries can also be used to create and initialize **mrsimulator** objects.
    To do this, the dictionary must use the object's attribute names as the key strings and be
    passed to a higher level object. Since a SpectralEvent object holds a list of
    TransitionQuery objects, the above code could have been written as

    .. plot::
        :context: close-figs

            # Both SymmetryQuery and TransitionQuery object as dict
            symm_query_dict = {"P": [-1], "D": [1]}
            trans_query_dict = {"ch1": symm_query_dict}

            # Dictionary of TransitionQuery passed to SpectralEvent
            spec_event = SpectralEvent(transition_queries=[trans_query_dict])

In the example above, the SymmetryQuery object is created and assigned to
the TransitionQuery attribute ``ch1``, i.e., it acts on the isotope in the
"first channel". Recall that the ``channels`` attribute of the Method object
holds an ordered list of isotope strings. This list's first, second, and third
isotopes are associated with ``ch1``, ``ch2``, and ``ch3``, respectively.
Currently, **mrsimulator** only supports up to three channels, although this may
be increased in future versions.

The TransitionQuery object goes into a list in the
``transition_queries`` attribute of a SpectralEvent object. The
SpectralEvent object, in turn, is added to an ordered list in the ``events``
attribute of a SpectralDimension object. All this is illustrated in the code
below.

.. plot::
    :context: close-figs

    from mrsimulator import Site, Coupling, SpinSystem, Simulator
    from mrsimulator import Method, SpectralDimension
    from mrsimulator import signal_processor as sp
    import matplotlib.pyplot as plt
    import numpy as np

    # Create single Site and Spin System
    deuterium = Site(
        isotope="2H",
        isotropic_chemical_shift=10,  # in ppm
        shielding_symmetric={"zeta": -80, "eta": 0.25},  # zeta in ppm
        quadrupolar={"Cq": 10e3, "eta": 0.0, "alpha": 0, "beta": np.pi / 2, "gamma": 0},
    )
    deuterium_system = SpinSystem(sites=[deuterium])

    # This method selects all observable (p_I = –1) transitions
    method_both_transitions = Method(
        channels=["2H"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            SpectralDimension(
                count=512,
                spectral_width=40000,  # in Hz
                events=[SpectralEvent(transition_queries=[{"ch1": {"P": [-1]}}])],
            )
        ],
    )

    # This method selects observable (p_I = –1) transitions with d_I = 1
    method_transition1 = Method(
        channels=["2H"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            SpectralDimension(
                count=512,
                spectral_width=40000,  # in Hz
                events=[SpectralEvent(transition_queries=[{"ch1": {"P": [-1], "D": [1]}}])],
            )
        ],
    )

    # This method selects observable (p_I = –1) transitions with d_I = -1
    method_transition2 = Method(
        channels=["2H"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            SpectralDimension(
                count=512,
                spectral_width=40000,  # in Hz
                events=[SpectralEvent(transition_queries=[{"ch1": {"P": [-1], "D": [-1]}}])],
            )
        ],
    )
    # Simulate spectra for all three method with spin system
    sim = Simulator(
        spin_systems=[deuterium_system],
        methods=[method_both_transitions, method_transition1, method_transition2],
    )
    sim.run()

    # Create SignalProcessor for Gaussian Convolution
    processor = sp.SignalProcessor(
        operations=[sp.IFFT(), sp.apodization.Gaussian(FWHM="100 Hz"), sp.FFT()]
    )

.. skip: next

.. plot::
    :context: close-figs

    # Plot spectra from all three methods
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.5), subplot_kw={"projection": "csdm"})
    ax[0].plot(
        processor.apply_operations(dataset=sim.methods[0].simulation).real,
        label="$p_I = -1$ transition",
    )
    ax[0].set_title("Single-Quantum Spectrum All Transitions")
    ax[0].legend()
    ax[0].grid()
    ax[0].invert_xaxis()  # reverse x-axis
    ax[1].plot(
        processor.apply_operations(dataset=sim.methods[1].simulation).real,
        label="$(p_I,d_I) = (-1,+1)$ transitions",
    )
    ax[1].plot(
        processor.apply_operations(dataset=sim.methods[2].simulation).real,
        label="$(p_I,d_I) = (-1,-1)$ transitions",
    )
    ax[1].set_title("Single-Quantum Spectrum Single Transitions")
    ax[1].legend()
    ax[1].grid()
    ax[1].invert_xaxis()  # reverse x-axis
    ax[0].set_ylim(-0.02, 0.34)  # Set y-limits to be the same
    ax[1].set_ylim(-0.02, 0.34)  # on both plots
    plt.tight_layout()
    plt.show()

.. note::

    Whenever the ``D`` attribute is omitted, the SymmetryQuery allows
    transitions with all values of :math:`\text{d}_I`. On the other hand, whenever
    the ``P`` attribute is omitted, it defaults to ``P=[0]``,  i.e., no selected
    transitions on the assigned channel.

Selecting Symmetric Single-Spin Transitions
'''''''''''''''''''''''''''''''''''''''''''

A notable case, particularly useful for half-integer quadrupolar nuclei, is that
:math:`\text{d}_I = 0` for all symmetric :math:`(m \rightarrow - m)`
transitions, as these transitions are unaffected by the first-order quadrupolar
coupling frequency contribution.  The MQ-MAS experiment involves a 2D correlation of
the two symmetric (:math:`\text{d}_I = 0`) transitions,
:math:`\ketbra{-\tfrac{1}{2}}{\tfrac{1}{2}}`, the so-called "central transition,"
and :math:`\ketbra{-\tfrac{3}{2}}{\tfrac{3}{2}}`, the symmetric triple quantum
transition.   The code below is an example of a custom 2D method using two
SpectralDimension objects, each holding a single SpectralEvent.  The
TransitionQuery objects select each transition in their respective
SpectralDimension objects.

.. plot::
    :context: close-figs

    my_mqmas = Method(
        channels=["87Rb"],
        magnetic_flux_density=9.4,
        rotor_frequency=np.inf,  # in Hz (here, set to infinity)
        spectral_dimensions=[
            SpectralDimension(
                count=128,
                spectral_width=6e3,  # in Hz
                reference_offset=-9e3,  # in Hz
                label="Symmetric 3Q Frequency",
                events=[SpectralEvent(transition_queries=[{"ch1": {"P": [-3], "D": [0]}}])],
            ),
            SpectralDimension(
                count=256,
                spectral_width=6e3,  # in Hz
                reference_offset=-5e3,  # in Hz
                label="Central Transition Frequency",
                events=[SpectralEvent(transition_queries=[{"ch1": {"P": [-1], "D": [0]}}])],
            ),
        ],
    )

    # Create three sites in RbNO3
    site1 = Site(
        isotope="87Rb",
        isotropic_chemical_shift=-27.4,  # ppm
        quadrupolar={"Cq": 1.68e6, "eta": 0.2},  # Cq in Hz
    )
    site2 = Site(
        isotope="87Rb",
        isotropic_chemical_shift=-28.5,  # ppm
        quadrupolar={"Cq": 1.94e6, "eta": 1},  # Cq in Hz
    )
    site3 = Site(
        isotope="87Rb",
        isotropic_chemical_shift=-31.3,  # ppm
        quadrupolar={"Cq": 1.72e6, "eta": 0.5},  # Cq in Hz
    )

    # No Couplings, so create a separate SpinSystem for each site.
    sites = [site1, site2, site3]
    RbNO3_spin_systems = [SpinSystem(sites=[s]) for s in sites]

    sim = Simulator(spin_systems=RbNO3_spin_systems, methods=[my_mqmas])
    sim.run()

    # Apply Gaussian line broadening along both dimensions via convolution
    gauss_convolve = sp.SignalProcessor(
        operations=[
            sp.IFFT(dim_index=(0, 1)),
            sp.apodization.Gaussian(FWHM="0.08 kHz", dim_index=0),
            sp.apodization.Gaussian(FWHM="0.22 kHz", dim_index=1),
            sp.FFT(dim_index=(0, 1)),
        ]
    )
    dataset = gauss_convolve.apply_operations(dataset=sim.methods[0].simulation)

In the code below, we use the `imshow()
<https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.imshow.html>`__ from the
matplotlib.pyplot module to
return an image of the dataset on a 2D regular raster. We also use ``"gist_ncar_r"`` from
`matpltolib's included colormaps
<https://matplotlib.org/stable/gallery/color/colormap_reference.html>`__
to map the dataset amplitude to colors; the `colorbar()
<https://matplotlib.org/stable/api/colorbar_api.html?highlight=colorbar#module-matplotlib.colorbar>`__
function provides the visualization of the dataset mapping to color to the right
of the plot.

.. skip: next

.. plot::
    :context: close-figs

    plt.figure(figsize=(4, 3))
    ax = plt.subplot(projection="csdm")
    cb = ax.imshow(dataset.real / dataset.real.max(), aspect="auto", cmap="gist_ncar_r")
    plt.colorbar(cb)
    ax.invert_xaxis()
    ax.invert_yaxis()
    plt.tight_layout()
    plt.show()

.. warning::

    This custom method, as well as the built-in Multi-Quantum VAS methods, assumes uniform
    excitation and mixing of the multiple-quantum transition. In an experimental MQ-MAS
    measurement, excitation and mixing efficiencies depend on the ratio of the quadrupolar
    coupling constant to the rf field strength. Therefore, the relative integrated intensities
    of this simulation may not agree with the experiment.


Inspecting Transition and Symmetry Pathways
'''''''''''''''''''''''''''''''''''''''''''

You can view the symmetry pathways that will be selected by your custom method
in a given spin system using the function
:py:meth:`~mrsimulator.Method.get_symmetry_pathways` as shown below.

.. plot::
    :context: close-figs

    from pprint import pprint
    pprint(my_mqmas.get_symmetry_pathways("P"))
    pprint(my_mqmas.get_symmetry_pathways("D"))

.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    [SymmetryPathway(
        ch1(87Rb): [-3] ⟶ [-1]
        total: -3.0 ⟶ -1.0
    )]
    [SymmetryPathway(
        ch1(87Rb): [0] ⟶ [0]
        total: 0.0 ⟶ 0.0
    )]


.. Method also has a related function :py:meth:`~mrsimulator.Method.plot` for generating a symmetry
.. pathway diagram of the method.
..
.. .. skip: next
..
.. .. plot::
..     :context: close-figs
..
..     pathway_diagram = my_mqmas.plot()
..     pathway_diagram.show()


Similarly, you can view the transition pathway that will be selected by your custom method in a given
spin system using the function :py:meth:`~mrsimulator.Method.get_transition_pathways` as shown
below.

.. plot::
    :context: close-figs

    from pprint import pprint
    pprint(my_mqmas.get_transition_pathways(SpinSystem(sites=[site1])))

.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    [|-1.5⟩⟨1.5| ⟶ |-0.5⟩⟨0.5|, weight=(1+0j)]


Multi-Spin Queries
------------------

When there is more than one site in a spin system, things get a little more
complicated with the SymmetryQuery objects.  Here we review some important
concepts associated with transition symmetry functions in coupled spin systems,
and see how SymmetryQuery objects are designed to work in such cases.

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
:math:`\hat{M}_2`, :math:`\hat{M}_3`, :math:`\hat{M}_4`, :math:`\hat{X}_1`,
:math:`\hat{X}_2`, :math:`\hat{X}_3`, and :math:`\hat{X}_4`
in the energy level diagram below.

.. figure:: ../../_static/ThreeCoupledSpinsEnergy.*
    :width: 550
    :alt: figure
    :align: center

    Energy level diagram for three coupled spin :math:`I=1/2` nuclei. Arrows
    beginning at the initial state and end at the final state represent the
    single-spin single-quantum transitions.   Transitions are labeled with their
    corresponding single-spin :math:`\text{p}_i` transition symmetry function
    values.

Keep in mind that the Method object does not know, in advance, the
number of sites in a spin system.

The TransitionQuery for selecting these 12 *single-spin single-quantum* transitions
is given in the code below.

.. plot::
    :context: close-figs

    # Create Site, Coupling and SpinSystem objects
    site_A = Site(isotope="1H", isotropic_chemical_shift=0.5)
    site_M = Site(isotope="1H", isotropic_chemical_shift=2.5)
    site_X = Site(isotope="1H", isotropic_chemical_shift=4.5)
    sites = [site_A, site_M, site_X]
    coupling_AM = Coupling(site_index=[0, 1], isotropic_j=12)
    coupling_AX = Coupling(site_index=[0, 2], isotropic_j=12)
    coupling_MX = Coupling(site_index=[1, 2], isotropic_j=12)
    couplings = [coupling_AM, coupling_AX, coupling_MX]
    proton_system = SpinSystem(sites=sites, couplings=couplings)

    # Custom method object emulating a one-pulse acquire experiment
    method = Method(
        channels=["1H"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            SpectralDimension(
                count=16000,
                spectral_width=1800,  # in Hz
                reference_offset=1000,  # in Hz
                label="$^{1}$H frequency",
                events=[SpectralEvent(transition_queries=[{"ch1": {"P": [-1]}}])],
            )
        ],
    )

    sim = Simulator(spin_systems=[proton_system], methods=[method])
    sim.run()

    # Add line broadening
    processor = sp.SignalProcessor(
        operations=[sp.IFFT(), sp.apodization.Exponential(FWHM="1 Hz"), sp.FFT()]
    )

.. skip: next

.. plot::
    :context: close-figs

    plt.figure(figsize=(10, 3))  # set the figure size
    ax = plt.subplot(projection="csdm")
    ax.plot(processor.apply_operations(dataset=sim.methods[0].simulation).real)
    ax.invert_xaxis()  # reverse x-axis
    plt.tight_layout()
    plt.grid()
    plt.show()

The assignment of transitions in the spectrum above are, from left to right, are
:math:`\hat{X}_4, (\hat{X}_3, \hat{X}_2)`, and :math:`\hat{X}_1` centered at 4.5
ppm, :math:`\hat{M}_4, (\hat{M}_3, \hat{M}_2)`, and :math:`\hat{M}_1` centered
at 2.5 ppm, and :math:`\hat{A}_4, (\hat{A}_3, \hat{A}_2)`, and :math:`\hat{A}_1`
centered at 0.5 ppm.

It is essential to realize that all sites having the same isotope are
"indistinguishable" to a TransitionQuery object. Recall that ``ch1`` is
associated with the first isotope in the list of isotope strings assigned to the
Method attribute ``channels``. When the TransitionQuery above is combined with
the SpinSystem object with three :math:`^1\text{H}` Sites, it must first expand
its SymmetryQuery into an intermediate set of spin-system-specific symmetry
queries, illustrated by each row in the table below.

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

The intermediate spin-system-specific symmetry query in each row selects a
subset of transitions from the complete set of transitions. The final set of
selected transitions is obtained from the union of transition subsets from each
spin-system-specific symmetry query.

The :py:meth:`~mrsimulator.Method.get_transition_pathways` function will allow
you to inspect the transitions selected by the Method
in terms of the initial and final Zeeman eigenstate quantum numbers.

.. plot::
    :context: close-figs

    from pprint import pprint
    pprint(method.get_transition_pathways(proton_system))

.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    [|-0.5, -0.5, -0.5⟩⟨-0.5, -0.5, 0.5|, weight=(1+0j),
    |-0.5, -0.5, -0.5⟩⟨-0.5, 0.5, -0.5|, weight=(1+0j),
    |-0.5, -0.5, 0.5⟩⟨-0.5, 0.5, 0.5|, weight=(1+0j),
    |-0.5, 0.5, -0.5⟩⟨-0.5, 0.5, 0.5|, weight=(1+0j),
    |-0.5, -0.5, -0.5⟩⟨0.5, -0.5, -0.5|, weight=(1+0j),
    |-0.5, -0.5, 0.5⟩⟨0.5, -0.5, 0.5|, weight=(1+0j),
    |0.5, -0.5, -0.5⟩⟨0.5, -0.5, 0.5|, weight=(1+0j),
    |-0.5, 0.5, -0.5⟩⟨0.5, 0.5, -0.5|, weight=(1+0j),
    |0.5, -0.5, -0.5⟩⟨0.5, 0.5, -0.5|, weight=(1+0j),
    |-0.5, 0.5, 0.5⟩⟨0.5, 0.5, 0.5|, weight=(1+0j),
    |0.5, -0.5, 0.5⟩⟨0.5, 0.5, 0.5|, weight=(1+0j),
    |0.5, 0.5, -0.5⟩⟨0.5, 0.5, 0.5|, weight=(1+0j)]


To further illustrate how the TransitionQuery and SymmetryQuery objects work in
a multi-site spin system, let's examine a few more examples in the case of three
weakly coupled proton sites.


Two-Spin Double-Quantum Transitions
'''''''''''''''''''''''''''''''''''

In this spin system there are six *two-spin double-quantum transitions* where
:math:`\text{p}_{AMX} = \text{p}_{A} + \text{p}_{M} + \text{p}_{X} = -2` and
another six *two-spin double-quantum transitions* where
:math:`\text{p}_{AMX} = \text{p}_{A} + \text{p}_{M} + \text{p}_{X} = +2`.  The
:math:`\text{p}_{AMX} = -2` transitions are illustrated in the energy-level diagram
below.

.. figure:: ../../_static/ThreeCoupledSpinsDoubleQuantum.*
    :width: 500
    :alt: figure
    :align: center

    Energy level diagram for three coupled spin :math:`I=1/2` nuclei. Arrows
    beginning at the initial state and end at the final state represent the
    two-spin double-quantum transitions.   Transitions are labeled with their
    corresponding single-spin :math:`\text{p}_i` transition symmetry function
    values.

The code below will select the six *two-spin double-quantum transitions* where
:math:`\text{p}_{AMX} = -2`.

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
                events=[SpectralEvent(transition_queries=[{"ch1": {"P": [-1, -1]}}])],
            )
        ],
    )

    sim = Simulator(spin_systems=[proton_system], methods=[method])
    sim.run()

.. skip: next

.. plot::
    :context: close-figs

    plt.figure(figsize=(10, 3))  # set the figure size
    ax = plt.subplot(projection="csdm")
    ax.plot(processor.apply_operations(dataset=sim.methods[0].simulation).real)
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
spin-system-specific symmetry queries illustrated in the table below.

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

Again, each row's intermediate spin-system-specific symmetry query is used to
select a subset of transitions from the complete set of transitions. The final
set of selected transitions is obtained from the union of transition subsets
from each spin-system-specific symmetry query.

.. plot::
    :context: close-figs

    from pprint import pprint
    pprint(method.get_transition_pathways(proton_system))

.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    [|-0.5, -0.5, -0.5⟩⟨-0.5, 0.5, 0.5|, weight=(1+0j),
    |-0.5, -0.5, -0.5⟩⟨0.5, -0.5, 0.5|, weight=(1+0j),
    |-0.5, -0.5, -0.5⟩⟨0.5, 0.5, -0.5|, weight=(1+0j),
    |-0.5, -0.5, 0.5⟩⟨0.5, 0.5, 0.5|, weight=(1+0j),
    |-0.5, 0.5, -0.5⟩⟨0.5, 0.5, 0.5|, weight=(1+0j),
    |0.5, -0.5, -0.5⟩⟨0.5, 0.5, 0.5|, weight=(1+0j)]


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
    :width: 500
    :alt: figure
    :align: center

    Energy level diagram for three coupled spin :math:`I=1/2` nuclei. Arrows
    beginning at the initial state and end at the final state represent the
    three spin single-quantum transitions.   Transitions are labeled with their
    corresponding single-spin :math:`\text{p}_i` transition symmetry function
    values.

The code below will select these *three-spin single-quantum transitions*.

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
                events=[SpectralEvent(transition_queries=[{"ch1": {"P": [-1, -1, +1]}}])]
            )
        ],
    )

    sim = Simulator(spin_systems=[proton_system], methods=[method])
    sim.run()

.. skip: next

.. plot::
    :context: close-figs

    plt.figure(figsize=(10, 3))  # set the figure size
    ax = plt.subplot(projection="csdm")
    ax.plot(processor.apply_operations(dataset=sim.methods[0].simulation).real)
    ax.invert_xaxis()  # reverse x-axis
    plt.tight_layout()
    plt.grid()
    plt.show()

The assignment of transitions in the spectrum above are, from left to right, are
:math:`\hat{S}_{3,AMX}`, :math:`\hat{S}_{2,AMX}`, and :math:`\hat{S}_{1,AMX}`

Again, combined with the three-site SpinSystem object, the SymmetryQuery is
expanded into the set of spin-system-specific symmetry queries illustrated in
the table below.


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

.. plot::
    :context: close-figs

    from pprint import pprint
    pprint(method.get_transition_pathways(proton_system))

.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    [|0.5, -0.5, -0.5⟩⟨-0.5, 0.5, 0.5|, weight=(1+0j),
    |-0.5, 0.5, -0.5⟩⟨0.5, -0.5, 0.5|, weight=(1+0j),
    |-0.5, -0.5, 0.5⟩⟨0.5, 0.5, -0.5|, weight=(1+0j)]

As you can surmise from the examples, the attributes of SymmetryQuery, ``P`` and
``D``, hold a list of single-spin transition symmetry function values, and the
length of the list is the desired number of spins that are involved in the
transition.

Heteronuclear multiple-spin transitions
'''''''''''''''''''''''''''''''''''''''

How does ``D`` fit into the multi-site SymmetryQuery story? Consider the
case of two coupled hydrogens, except we replace one of the :math:`^1H` with
:math:`^2H`.  Let's focus on the single-spin single-quantum transitions, shown
below as :math:`\hat{A}_{1\pm}` and :math:`\hat{A}_{2\pm}` on the left, and the
two-spin triple-quantum transition, shown below as  :math:`\hat{T}_{AX}` on
the right.

.. figure:: ../../_static/Spin1SpinHalfCouple.*
    :width: 850
    :alt: figure
    :align: center

    Energy level diagram for two coupled nuclei with spins :math:`I=1/2` and
    :math:`I=1`. Arrows beginning at the initial state and end at the final
    state represent the single spin single-quantum transitions (left) and the
    three-spin triple-quantum transition.  Transitions are labeled with their
    corresponding single-spin :math:`\text{p}_i` transition symmetry function
    values.

.. plot::
    :context: close-figs

    from mrsimulator.spin_system.tensors import SymmetricTensor
    import numpy as np

    site_A = Site(
        isotope="2H",
        isotropic_chemical_shift=0.5,
        quadrupolar=SymmetricTensor(
            Cq=100000,  # in Hz
            eta=0.2,
            alpha=5 * np.pi / 180,
            beta=np.pi / 2,
            gamma=70 * np.pi / 180,
        ),
    )
    site_X = Site(isotope="1H", isotropic_chemical_shift=4.5)
    sites = [site_A, site_X]
    coupling_AX = Coupling(site_index=[0, 1], dipolar={"D": -20000})
    couplings = [coupling_AX]
    system_AX = SpinSystem(sites=sites, couplings=couplings)

    methodAll1Q = Method(
        channels=["2H", "1H"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            SpectralDimension(
                count=16000,
                spectral_width=200000,  # in Hz
                reference_offset=0,  # in Hz
                label="$^{2}$H frequency",
                events=[SpectralEvent(transition_queries=[{"ch1": {"P": [-1]}}])],
            )
        ],
    )

    methodHalf1Q = Method(
        channels=["2H", "1H"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            SpectralDimension(
                count=16000,
                spectral_width=200000,  # in Hz
                reference_offset=0,  # in Hz
                label="$^{2}$H frequency",
                events=[SpectralEvent(transition_queries=[{"ch1": {"P": [-1], "D": [-1]}}])],
            )
        ],
    )

    method3Q = Method(
        channels=["2H", "1H"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            SpectralDimension(
                count=16000,
                spectral_width=10000,  # in Hz
                reference_offset=5000,  # in Hz
                label="$^{2}$H frequency",
                events=[
                    SpectralEvent(
                        transition_queries=[{"ch1": {"P": [-2]}, "ch2": {"P": [-1]}}]
                    )
                ],
            )
        ],
    )
    processor = sp.SignalProcessor(
        operations=[sp.IFFT(), sp.apodization.Gaussian(FWHM="100 Hz"), sp.FFT()]
    )

    sim = Simulator(spin_systems=[system_AX], methods=[methodAll1Q, methodHalf1Q, method3Q])
    sim.config.integration_volume = "hemisphere"
    sim.run()

.. skip: next

.. plot::
    :context: close-figs

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.5), subplot_kw={"projection": "csdm"})
    ax[0].plot(processor.apply_operations(dataset=sim.methods[0].simulation).real)
    ax[0].set_title("Full Single-Quantum Spectrum")
    ax[0].grid()
    ax[0].invert_xaxis()  # reverse x-axis
    ax[1].plot(processor.apply_operations(dataset=sim.methods[1].simulation).real)
    ax[1].set_title("Half Single-Quantum Spectrum")
    ax[1].grid()
    ax[1].invert_xaxis()  # reverse x-axis
    plt.tight_layout()
    plt.show()

The deuterium spectrum of a static-polycrystalline sample is shown on the left
for all single-spin single-quantum transitions on deuterium,
:math:`\hat{A}_{1\pm}` and :math:`\hat{A}_{2\pm}`. The spectrum on the right is
for half of the single-spin single-quantum transitions on deuterium:
:math:`\hat{A}_{1-}` and :math:`\hat{A}_{2-}`.

.. skip: next

.. plot::
    :context: close-figs

    plt.figure(figsize=(10, 3))  # set the figure size
    ax = plt.subplot(projection="csdm")
    ax.set_title("Heteronuclear Two-Spin ($^2$H-$^1$H) Triple-Quantum Spectrum")
    ax.plot(processor.apply_operations(dataset=sim.methods[2].simulation).real)
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

The single transition in the heteronuclear two-spin
(:math:`^2\text{H}`-:math:`^1\text{H}`) triple-quantum spectrum is unaffected by
the dipolar and quadrupolar frequency anisotropies.


Frequency Contributions
-----------------------

The NMR frequency, :math:`\Omega(\Theta,i,j)`, of an :math:`i  \rightarrow  j`
transition between the eigenstates of the stationary-state semi-classical
Hamiltonian in a sample with a lattice spatial orientation, :math:`\Theta`, can
be written as a sum of components,

.. math::
    \Omega(\Theta,i,j) = \sum_k \Omega_k(\Theta,i,j)

with each component, :math:`\Omega_k(\Theta,i,j)`, separated into three parts:

.. math::
    \Omega_k(\Theta,i,j) = \omega_k \, {\Xi}^{(k)}_L (\Theta) \,{\xi}^{(k)}_\ell (i,j),

where :math:`{\xi}^{(k)}_\ell(i,j)` are the spin transition symmetry functions
described earlier, :math:`{\Xi}^{(k)}_L(\Theta)` are the spatial symmetry
functions, and :math:`\omega_k` gives the size of the kth frequency component.
The experimentalist indirectly influences a frequency component :math:`\Omega_k`
by direct manipulation of the quantum transition, :math:`i \rightarrow  j`, and
the spatial orientation,  :math:`\Theta` of the sample.

The function symbol :math:`\Xi_\ell(\Theta)` is replaced with the
upper-case symbols :math:`\mathbb{S}`, :math:`\mathbb{P}(\Theta)`,
:math:`\mathbb{D}(\Theta)`, :math:`\mathbb{F}(\Theta)`,
:math:`\mathbb{G}(\Theta)`, :math:`\ldots`, i.e., following the spectroscopic
sub-shell letter designations for :math:`L`. Consult the `Symmetry Pathways
paper <https://doi.org/10.1016/j.pnmrs.2010.11.003>`_ for more details on the
form of the spatial symmetry functions.  In short, the :math:`\mathbb{S}`
function is independent of sample orientation, i.e., it will appear in all
isotropic frequency contributions.  The :math:`\mathbb{D}(\Theta)` function has
a second-rank dependence on sample orientation, and can be averaged away with
fast magic-angle spinning, i.e., spinning about an angle, :math:`\theta_R`, that
is the root of the second-rank Legendre polynomial :math:`P_2(\cos \theta_R)`.
The other spatial symmetry functions are removed by spinning the sample about
the corresponding root of the :math:`L`th-rank Legendre polynomial ":math:`P_L(\cos
\theta_R)`.

.. note::

    For 2nd-order quadrupolar coupling contributions, it is convenient to define
    "hybrid" spin transition functions as linear combinations of the spin transition
    functions

    .. math::

        \mathbb{c}_0  = \,\,\,\frac{4}{\sqrt{125}} \, [I(I+1) - 3/4] \, \mathbb{p}_I  + \sqrt{\frac{18}{25}} \, \mathbb{f}_I

    .. math::

        \mathbb{c}_2  = \,\,\,\frac{2}{\sqrt{175}} \, [I(I+1) - 3/4] \, \mathbb{p}_I  - \frac{6}{\sqrt{35}} \, \mathbb{f}_I

    .. math::

        \mathbb{c}_4  = -\frac{184}{\sqrt{875}} \, [I(I+1) - 3/4] \, \mathbb{p}_I  - \frac{17}{\sqrt{175}} \, \mathbb{f}_I

These transition symmetry functions play an essential role in evaluating the
individual frequency contributions to the overall transition frequency, given in
the table below and in
:py:meth:`~mrsimulator.method.frequency_contrib.FrequencyEnum`. They also aid in
pulse sequence design by identifying how different frequency contributions
refocus through the transition pathways.  For a summary on echo symmetry classification in NMR,
click on the disclosure button below.

.. raw:: html

    <details class="reveal-theory">
    <summary><b>Echo Symmetry Classification</b></summary>

.. include:: echo_symm_classes.rst

.. raw:: html

    </details>


.. _frequency_contribution_table:

.. list-table:: Frequency Contributions
    :widths: 25 25 25 25 25
    :header-rows: 2

    * - Interactions
      - perturbation
      - anisotropy
      - ``freq_contrib``
      - Expression
    * -
      - order
      - rank
      -
      -
    * - shielding
      - 1st
      - 0th
      - ``Shielding1_0``
      - :math:`-\omega_0 \sigma_\text{iso} \cdot \mathbb{p}_I`
    * - shielding
      - 1st
      - 2nd
      - ``Shielding1_2``
      - :math:`-\omega_0 \zeta_\sigma \cdot \mathbb{D}^{\{\sigma\}} \cdot \mathbb{p}_I`
    * - weak J
      - 1st
      - 0th
      - ``J1_0``
      - :math:`2 \pi J_\text{iso} \, (\mathbb{pp})_{IS}`
    * - weak J
      - 1st
      - 2nd
      - ``J1_2``
      - :math:`2 \pi \zeta_J \cdot \mathbb{D}^{\{d_{IS}\}} \cdot (\mathbb{pp})_{IS}`
    * - weak dipolar
      - 1st
      - 2nd
      - ``D1_2``
      - :math:`\omega_d \cdot \mathbb{D}^{\{d_{IS}\}} \cdot (\mathbb{pp})_{IS}`
    * - quadrupolar
      - 1st
      - 2nd
      - ``Quad1_2``
      - :math:`\omega_q \cdot \mathbb{D}^{\{q\}} \cdot \mathbb{d}_I`
    * - quadrupolar
      - 2nd
      - 0th
      - ``Quad2_0``
      - :math:`\displaystyle \frac{\omega_q^2}{\omega_0}  \cdot \mathbb{S}^{\{qq\}} \cdot \mathbb{c}_0`
    * - quadrupolar
      - 2nd
      - 2nd
      - ``Quad2_2``
      - :math:`\displaystyle\frac{\omega_q^2}{\omega_0}  \cdot \mathbb{D}^{\{qq\}} \cdot \mathbb{c}_2`
    * - quadrupolar
      - 2nd
      - 4th
      - ``Quad2_4``
      - :math:`\displaystyle\frac{\omega_q^2}{\omega_0}  \cdot \mathbb{G}^{\{qq\}} \cdot \mathbb{c}_4`


Affine Transformations
----------------------

The ability to refocus different spatial and transition symmetries into echoes
with different paths in time-resolved NMR experiments creates opportunities for
generating multi-dimensional spectra that correlate different interactions.
These spectra can be made easier to interpret through similarity
transformations. Most similarity transformations in NMR are affine
transformations, as they preserve the colinearity of points and ratios of
distances. Essential in any similarity transformation is whether to implement
the transformation actively or passively. Active transformations change the
appearance of the signal while leaving the coordinate system unchanged, whereas
passive transformations leave the appearance of the signal unchanged while
changing the coordinate system. Both active and passive transformations are used
extensively in NMR.

The general form of the affine transformation of a n-dimensional spectrum is

.. math::

    {\boldsymbol \Omega}' = {\cal A} {\boldsymbol \Omega}

In the two-dimensional case, this is given by

.. math::
    \left[
    \begin{array}{c}
    \Omega^{'[1]} \\
    \Omega^{'[2]}
    \end{array}
    \right]
    =
    \underbrace{
    \left[
    \begin{array}{cc}
    a & b \\
    c & d
    \end{array}
    \right]
    }_{\cal A}
    \left[
    \begin{array}{c}
    \Omega^{[1]} \\
    \Omega^{[2]}
    \end{array}
    \right]

.. note::

    For the multiple-quantum MAS experiment, a shear and scale transformation is
    often applied to the spectrum to create a 2D spectrum correlating the MQ-MAS
    isotropic frequency to the anisotropic central transition frequency. This
    correlation can be achieved by adding an affine matrix to the method.

    For 3Q-MAS on a spin :math:`I=3/2` nucleus, where the shear factor is
    :math:`\kappa^{(\omega_2)} = 21/27`, the affine matrix giving the
    appropriate shear and scale transformation is given by

    .. math::
        {\cal A}_2 =
        \left[
        \begin{array}{cc}
        \displaystyle \frac{1}{1 + |\kappa^{(\omega_2)}|}
        & \displaystyle \frac{	\kappa^{(\omega_2)}}{1 + |\kappa^{(\omega_2)}| } \\
        0 & 1
        \end{array}
        \right]
        =
        \left[
        \begin{array}{cc}
        9/16 & 7/16 \\
        0 & 1
        \end{array}
        \right]

    After the affine transformation, the position of the resonance in the
    isotropic projection is a weighted average of the multiple quantum and
    central transition isotropic frequencies given by

    .. math::
        \left \langle\Omega_{iso} \right \rangle_{\text{MQ-MAS}}
        =
        \frac{1}{1 + |\kappa^{(\omega_1)}|}
        \,
        \Omega_\text{iso}(m,-m)
        +
        \frac{\kappa^{(\omega_1)}}{1 + |\kappa^{(\omega_1)}|}
        \,
        \Omega_\text{iso}\left(\textstyle \frac{1}{2},-\frac{1}{2}\right).

    If the spectrum is to be referenced to a frequency other than the rf carrier
    frequency (i.e. zero is not defined in the middle of the spectrum), then the
    reference offset used in the single-quantum dimension must be multiplied by a
    factor of
    :math:`{\left({\text{p}_I^{[1]}}/{\text{p}_I^{[2]}} + |\kappa^{(\omega_1)}| \right)/(1+ |\kappa^{(\omega_1)}| )}`
    when used in the isotropic dimension.

    See the `"Symmetry Pathways in Solid-State NMR" paper
    <https://doi.org/10.1016/j.pnmrs.2010.11.003>`_  for a more detailed
    discussion on affine transformations in NMR.

In the code below, the 3Q-MAS method described earlier is modified to include an
affine matrix to perform this shear transformation.

.. plot::
    :context: close-figs

    my_sheared_mqmas = Method(
        channels=["87Rb"],
        magnetic_flux_density=9.4,
        rotor_frequency=np.inf,  # in Hz (here, set to infinity)
        spectral_dimensions=[
            SpectralDimension(
                count=128,
                spectral_width=6e3,  # in Hz
                reference_offset=-9e3,  # in Hz
                label="3Q-MAS isotropic dimension",
                events=[
                    SpectralEvent(transition_queries=[{"ch1": {"P": [-3], "D": [0]}}])
                ],
            ),
            SpectralDimension(
                count=256,
                spectral_width=6e3,  # in Hz
                reference_offset=-5e3,  # in Hz
                label="Central Transition Frequency",
                events=[
                    SpectralEvent(transition_queries=[{"ch1": {"P": [-1], "D": [0]}}])
                ],
            ),
        ],
        affine_matrix=[[9 / 16, 7 / 16], [0, 1]],
    )

    sim = Simulator(spin_systems=RbNO3_spin_systems, methods=[my_sheared_mqmas])
    sim.run()

    dataset = gauss_convolve.apply_operations(dataset=sim.methods[0].simulation)

.. skip: next

.. plot::
    :context: close-figs

    plt.figure(figsize=(4, 3))
    ax = plt.subplot(projection="csdm")
    cb = ax.imshow(dataset.real / dataset.real.max(), aspect="auto", cmap="gist_ncar_r")
    plt.colorbar(cb)
    ax.invert_xaxis()
    ax.invert_yaxis()
    plt.tight_layout()
    plt.show()

.. note::

    For MQ-MAS, a second shear and scale can be applied to remove isotropic
    chemical shift component along the :math:`\Omega^{[2]''}` axis.  For a
    spin :math:`I=3/2` nucleus, with a second shear factor of
    :math:`\kappa^{(\omega_1)} = - 8/17`, the affine matrix is given by

    .. math::
        {\cal A}_1 =
        \left[
        \begin{array}{cc}
        1 & 0 \\
        \displaystyle \frac{	\kappa^{(\omega_1)}}{1 + |\kappa^{(\omega_1)}| }
        & \displaystyle \frac{1}{1 + |\kappa^{(\omega_1)}|}
        \end{array}
        \right]
        =
        \left[
        \begin{array}{cc}
        1 & 0 \\
        -8/25 & 17/25
        \end{array}
        \right],

    and the product of the two affine transformations is

    .. math::
        {\cal A}_T = {\cal A}_1 {\cal A}_2
        =
        \left[
        \begin{array}{cc}
        1 & 0 \\
        -8/25 & 17/25
        \end{array}
        \right]
        \left[
        \begin{array}{cc}
        9/16 & 7/16 \\
        0 & 1
        \end{array}
        \right]
        =
        \left[
        \begin{array}{cc}
        9/16 & 7/16 \\
        -9/50 & 27/50
        \end{array}
        \right].

Below is the code for simulating a 3Q-MAS spectrum with a double shear transformation.

.. plot::
    :context: close-figs

    my_twice_sheared_mqmas = Method(
        channels=["87Rb"],
        magnetic_flux_density=9.4,
        rotor_frequency=np.inf,  # in Hz (here, set to infinity)
        spectral_dimensions=[
            SpectralDimension(
                count=128,
                spectral_width=6e3,  # in Hz
                reference_offset=-9e3,  # in Hz
                label="3Q-MAS isotropic dimension",
                events=[
                    SpectralEvent(transition_queries=[{"ch1": {"P": [-3], "D": [0]}}])
                ],
            ),
            SpectralDimension(
                count=256,
                spectral_width=6e3,  # in Hz
                reference_offset=0,  # in Hz
                label="CT Quad-Only Frequency",
                events=[
                    SpectralEvent(transition_queries=[{"ch1": {"P": [-1], "D": [0]}}])
                ],
            ),
        ],
        affine_matrix=[[9 / 16, 7 / 16], [-9 / 50, 27 / 50]],
    )

    sim = Simulator(spin_systems=RbNO3_spin_systems, methods=[my_twice_sheared_mqmas])
    sim.run()

    dataset = gauss_convolve.apply_operations(dataset=sim.methods[0].simulation)

.. skip: next

.. plot::
    :context: close-figs

    plt.figure(figsize=(4, 3))
    ax = plt.subplot(projection="csdm")
    cb = ax.imshow(dataset.real / dataset.real.max(), aspect="auto", cmap="gist_ncar_r")
    plt.colorbar(cb)
    ax.invert_xaxis()
    ax.invert_yaxis()
    plt.tight_layout()
    plt.show()


Average Frequency & Multiple Events
-----------------------------------

To illustrate the versatility of the Method object, we can also design an MQ-MAS
method that correlates the isotropic MQ-MAS frequency to the central transition
without the need for an affine transformation.  Recall that the 3Q-MAS isotropic
frequency on spin :math:`I=3/2` is given by

.. math::
    \Omega_\text{iso} =  \frac{9}{16}\Omega_{3Q} + \frac{7}{16}\Omega_{CT}.

As we saw at the beginning of this section, the first spectral dimension
derives its *average frequency*, :math:`\overline{\Omega}_1`, from a weighted
average of multiple transition frequencies. Thus, this weighted average frequency can
be obtained through the use of multiple SpectralEvent objects in  the
SpectralDimension associated with the isotropic dimension, as shown in the
code below.

.. skip: next

.. plot::
    :context: close-figs

    my_three_event_mqmas = Method(
        channels=["87Rb"],
        magnetic_flux_density=9.4,
        rotor_frequency=np.inf,  # in Hz (here, set to infinity)
        spectral_dimensions=[
            SpectralDimension(
                count=128,
                spectral_width=6e3,  # in Hz
                reference_offset=-9e3,  # in Hz
                label="3Q-MAS isotropic dimension",
                events=[
                    SpectralEvent(
                        fraction=9 / 16, transition_queries=[{"ch1": {"P": [-3], "D": [0]}}]
                    ),
                    SpectralEvent(
                        fraction=7 / 16, transition_queries=[{"ch1": {"P": [-1], "D": [0]}}]
                    ),
                ],
            ),
            SpectralDimension(
                count=256,
                spectral_width=6e3,  # in Hz
                reference_offset=-5e3,  # in Hz
                label="Central Transition Frequency",
                events=[
                    SpectralEvent(transition_queries=[{"ch1": {"P": [-1], "D": [0]}}])
                ],
            ),
        ],
    )

    sim = Simulator(spin_systems=RbNO3_spin_systems, methods=[my_three_event_mqmas])
    sim.run()

    dataset = gauss_convolve.apply_operations(dataset=sim.methods[0].simulation)

.. skip: next

.. plot::
    :context: close-figs

    plt.figure(figsize=(4, 3))
    ax = plt.subplot(projection="csdm")
    cb = ax.imshow(dataset.real / dataset.real.max(), aspect="auto", cmap="gist_ncar_r")
    plt.colorbar(cb)
    ax.invert_xaxis()
    ax.invert_yaxis()
    plt.tight_layout()
    plt.show()

We could apply an affine transformation to remove the isotropic chemical shift
from the central transition (horizontal) dimension. If you go back to the
previous discussion, you will find that the required value for the
``affine_matrix`` in the Method object to do this shear is given by

``affine_matrix=[[1,0],[-8/25, 17/25]]``

Mixing Queries
--------------

The amplitude of a transition pathway signal derives from the product
of mixing amplitudes associated with each transfer between transitions in a
transition pathway,
e.g.,

.. math::
    (a_{0A})\,\hat{t}_A \rightarrow (a_{0A}a_{AB})\,\hat{t}_B\rightarrow (a_{0A}a_{AB}a_{BC})\,\hat{t}_C \rightarrow \cdots

Here, :math:`a_{0A}` is the amplitude of the initial :math:`\hat{t}_A` transition,
:math:`a_{AB}` is the mixing amplitude for the transfer from
:math:`\hat{t}_A` to :math:`\hat{t}_B`,  :math:`a_{BC}` is the mixing amplitude
for the transfer from :math:`\hat{t}_B` to :math:`\hat{t}_C`, and so on.  The
growing product :math:`(a_{0A}a_{AB}a_{BC} \cdots)` is the transition pathway
amplitude.  Eliminating a transition with a TransitionQuery in a spectral or
delay event sets the eliminated transition's pathway amplitude to zero, i.e., it
prunes that transition pathway branch.

Default Total Mixing between Adjacent Spectral or Delay Events
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

In previous discussions, we did not mention the efficiency of transfer between
selected transitions in adjacent SpectralEvent objects. This is because, as
default behavior, **mrsimulator** does a *total mixing*, i.e., connects all
selected transitions in the two adjacent spectral or delay events. In other
words, if the first of two adjacent SpectralEvent objects has three selected
transitions, and the second has two selected transitions, then **mrsimulator**
will make :math:`3 \times 2 = 6` connections, i.e., six transition pathways
passing from the first to second SpectralEvent objects.

Additionally, this *total mixing* assumes that every connection has a mixing
amplitude of 1. This is unrealistic, but if used correctly gives a significant
speed-up in the simulation by avoiding the need to calculate mixing amplitudes.

.. warning::

    The use of total mixing, i.e., the default mixing, can complicate the comparison
    of integrated intensities between different methods, depending on the selected
    transition pathways.

If this default mixing behavior had been explicitly shown in the previous example,
the events list in the first SpectralDimension would have looked like the code
below.

.. plot::
    :context: close-figs

    from mrsimulator.method import MixingEvent

    events = [
        SpectralEvent(fraction=9 / 16, transition_queries=[{"ch1": {"P": [-3], "D": [0]}}]),
        MixingEvent(query="TotalMixing"),
        SpectralEvent(fraction=7 / 16, transition_queries=[{"ch1": {"P": [-1], "D": [0]}}]),
        MixingEvent(query="TotalMixing"),
    ]

Since only one transition was selected in each SpectralEvent, the expected (and
default) behavior is that there is a mixing (transfer) of coherence between
the symmetric triple-quantum and central transitions, forming the desired
transition pathway.

However, when multiple transition pathways are present in a method, you may need
more accurate mixing amplitudes when connecting selected transitions of adjacent
events. You may also need to prevent the undesired mixing of specific
transitions between two adjacent events. As described below, you can avoid a
``"TotalMixing"`` event by inserting MixingEvent object with a certain rotation
query.

Rotation Query
''''''''''''''

A rotation of :math:`\theta` about an axis defined by :math:`\phi`  in the
:math:`x`-:math:`y` plane on a selected transition, :math:`\ketbra{I, m_f}{I,
m_i}`, in a spectral or delay event transfers it to all selected transitions,
:math:`\ketbra{I,m_f'}{I,m_i'}` in the next spectral or delay event, according
to

.. math::

    \ketbra{I, m_f}{I, m_i} \stackrel{\theta_\phi}{\longrightarrow} \sum_{m_f'}\sum_{m_i'} d_{m_f',m_f}^{(I)}(\theta)d_{m_i',m_i}^{(I)}(\theta)e^{-i\Delta p\phi}(i)^{\Delta p}\ketbra{I,m_f'}{I,m_i'},

where :math:`\Delta p_I = p_I' - p_I`.  From this expression, we obtain the complex mixing amplitude
from :math:`\ketbra{I, m_f}{I, m_i}` to :math:`\ketbra{I, m_f'}{I, m_i'}` due to a rotation to be

.. math::

    a(\theta,\phi) = d_{m_f',m_f}^{(I)}(\theta)d_{m_i',m_i}^{(I)}(\theta)e^{-i\Delta p\phi}(i)^{\Delta p}.

From this expression, we note a few interesting and useful cases.  One is the coherence
transfer under a :math:`\pi` rotation, given by

.. math::
    :label: piPulseTransition

    \ketbra{I,m_f}{I, m_i}  \stackrel{\pi_\phi}{\longrightarrow} \ketbra{I, -m_f}{I, -m_i} e^{-i\Delta p\phi}(i)^{\Delta p}.

That is, a :math:`\pi` rotation will make only one connection between
transitions in adjacent events.  It is also a special connection because the
:math:`\text{p}_I` transition symmetry value for the two transitions are equal
but opposite in sign.  Additionally, the :math:`\text{d}_I` transition symmetry
remains unchanged (:math:`\Delta \text{d}_I = 0`) for the two transitions. (In
fact, this behavior under a :math:`\pi` rotation is generally true for odd
(:math:`\text{p}_I, \text{f}_I, \ldots)` and even (:math:`\text{d}_I,
\text{g}_I, \ldots)` rank spin transition symmetry functions.)

Another interesting result is that, while a rotation can transfer a transition
into many other transitions, the :math:`\text{d}_I` transition symmetry value
cannot remain unchanged (:math:`\Delta \text{d}_I \neq 0`) between two connected
transitions under a :math:`\pi/2` rotation.

Finally, another useful result is

.. math::
    :label: zeroPulseTransition

    \ketbra{I,m_f}{I, m_i}  \stackrel{(0)_\phi}{\longrightarrow} \ketbra{I, m_f}{I, m_i}.

While it's not surprising that a rotation through an angle of zero does nothing
to the transition, this turns out to help act as the opposite of a total mixing
event, i.e., a ``"NoMixing"`` event. As a convenience, this is defined as a
``"NoMixing"`` query and can be implemented with the code below.

.. plot::
    :context: close-figs

    MixingEvent(query="NoMixing")

The MixingEvent object holds the rotation details in a MixingQuery object as
a RotationQuery object associated with a ``channels`` attribute.  This is
illustrated in the sample code below.

.. plot::
    :context: close-figs

    import numpy as np
    from mrsimulator.method.query import RotationQuery
    rot_query_90 = RotationQuery(angle=np.pi/2, phase=0)
    rot_query_180 = RotationQuery(angle=np.pi, phase=0)
    rot_mixing = MixingEvent(query={
            "ch1": rot_query_90,
            "ch2": rot_query_180
        }
    )


p and d Echoes on Deuterium
'''''''''''''''''''''''''''

Here, we examine two examples in a deuterium spin system that illustrate the
importance of echo classification in understanding how transition-frequency
contributions can be eliminated or separated based on their dependence on
different transition symmetry functions.

First, we implement two Method objects that follow the design of the
experimental pulse sequence. In this effort, we use RotationQuery objects to
select the desired transition pathways and obtain spectra with the desired
average frequencies. Then, we implement two simpler Method objects that
produce identical spectra and illustrate how :ref:`frequency
contributions<freq_contrib_api>` can be used to reduce the number of events
needed in a custom method.

Consider the Hahn and Solid Echo pulse sequences on the left and right,
respectively.

.. figure:: ../../_static/HahnAndSolidEcho.*
    :alt: Transition symmetry pathways for the Hahn and Solid Echo experiments
    :align: center
    :width: 100%

    Hahn Echo (left) and Solid-Echo (right) pulse sequences. Above each
    sequence, on the energy level diagram, are the corresponding two transition
    pathways indicated with blue and orange arrows.  Transitions are labeled
    with their corresponding :math:`\text{p}_I` and :math:`\text{d}_I`
    transition symmetry function values.  Below each sequence are the
    corresponding transition symmetry pathways, also in blue and orange.

The Hahn Echo sequence, with :math:`\pi/2-\tau-\pi-t\rightarrow`, leads to the formation
of a :math:`\text{p}_I` echo at :math:`t = \tau`.  The two transition pathways
created by this experiment on a deuterium nucleus are illustrated beneath the
sequence. Remember that a :math:`\pi` rotation is a special because it connects
transitions with equal but opposite signs of :math:`\text{p}_I` while
:math:`\text{d}_I` remains invariant.

The Solid Echo sequence, with :math:`\pi/2-\tau-\pi/2-t\rightarrow`, leads to the
formation of a :math:`\text{d}_I` echo at :math:`t = \tau`.  The two transition
pathways created by this experiment on a deuterium nucleus are illustrated
beneath the sequence. Here, also recall that the :math:`\text{d}_I` transition
symmetry value cannot remain unchanged (:math:`\Delta \text{d}_I \neq 0`)
between two connected transitions under a :math:`\pi/2` rotation.

Below are two custom Method objects for simulating the Hahn and Solid Echo
experiments. There is only one SpectralDimension object in each method, and
the average frequency during each spectral dimension is derived from equal
fractions of two SpectralEvent objects.  Between these two SpectralEvent
objects is a MixingEvent with a RotationQuery object. The
RotationQuery object is created with a :math:`\pi` rotation in the Hahn Echo
method, and a :math:`\pi/2` rotation in the Solid Echo method.

.. note ::

    The ``transition_queries`` attribute of SpectralEvent holds a list of
    TransitionQuery objects. Each TransitionQuery in the list applies to
    the full set of transitions in the spin system. The union of these transition
    subsets becomes the final set of selected transitions during the
    SpectralEvent.

We use the deuterium Site defined earlier in this document.

.. plot::
    :context: close-figs

    from mrsimulator.method import MixingEvent

    deuterium = Site(
        isotope="2H",
        isotropic_chemical_shift=10,  # in ppm
        shielding_symmetric={"zeta": -80, "eta": 0.25},  # zeta in ppm
        quadrupolar={"Cq": 10e3, "eta": 0.0, "alpha": 0, "beta": np.pi / 2, "gamma": 0},
    )
    deuterium_system = SpinSystem(sites=[deuterium])

    hahn_echo = Method(
        channels=["2H"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            SpectralDimension(
                count=512,
                spectral_width=2e4,  # in Hz
                events=[
                    SpectralEvent(
                        fraction=0.5,
                        transition_queries=[
                            {"ch1": {"P": [1], "D": [1]}},
                            {"ch1": {"P": [1], "D": [-1]}},
                        ],
                    ),
                    MixingEvent(query={"ch1": {"angle": 3.141592, "phase": 0}}),
                    SpectralEvent(
                        fraction=0.5,
                        transition_queries=[
                            {"ch1": {"P": [-1], "D": [1]}},
                            {"ch1": {"P": [-1], "D": [-1]}},
                        ],
                    ),
                ],
            )
        ],
    )

    solid_echo = Method(
        channels=["2H"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            SpectralDimension(
                count=512,
                spectral_width=2e4,  # in Hz
                events=[
                    SpectralEvent(
                        fraction=0.5,
                        transition_queries=[
                            {"ch1": {"P": [-1], "D": [1]}},
                            {"ch1": {"P": [-1], "D": [-1]}},
                        ],
                    ),
                    MixingEvent(query={"ch1": {"angle": 3.141592 / 2, "phase": 0}}),
                    SpectralEvent(
                        fraction=0.5,
                        transition_queries=[
                            {"ch1": {"P": [-1], "D": [1]}},
                            {"ch1": {"P": [-1], "D": [-1]}},
                        ],
                    ),
                ],
            )
        ],
    )

We can check the resulting transition pathways using these TransitionQuery objects with the
code below for the ``hahn_echo`` method,

.. plot::
    :context: close-figs

    from pprint import pprint
    pprint(hahn_echo.get_transition_pathways(deuterium_system))

.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    [|1.0⟩⟨0.0| ⟶ |-1.0⟩⟨0.0|, weight=(1+0j)
     |0.0⟩⟨-1.0| ⟶ |0.0⟩⟨1.0|, weight=(1+0j)]

and for the ``solid_echo`` method with the code below.

.. plot::
    :context: close-figs

    pprint(solid_echo.get_transition_pathways(deuterium_system))

.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    [|-1.0⟩⟨0.0| ⟶ |0.0⟩⟨1.0|, weight=(0.5+0j)
     |0.0⟩⟨1.0| ⟶ |-1.0⟩⟨0.0|, weight=(0.5+0j)]

Notice that the weights of the transition pathways in the solid-echo method are
half of those in the Hahn-echo method. This is because the :math:`\pi` pulse in
the Hahn-echo method gives perfect transfer between the two transitions in the
adjacent spectral events. In contrast, while the :math:`\pi/2` pulse in the
solid-echo method prevents the undesired transition pathways with :math:`\Delta
\text{d}_I = 0`, it also connects the selected transitions during the first
spectral event to undesired transitions in the second spectral event, which are
eliminated by its symmetry query.

Next, we simulate both methods, and perform a Gaussian line shape convolution on
each output spectrum, and plot the datasets.

.. plot::
    :context: close-figs

    sim = Simulator()
    sim.spin_systems = [deuterium_system]
    sim.methods = [hahn_echo, solid_echo]
    sim.run()

    processor = sp.SignalProcessor(
        operations=[
            sp.IFFT(),
            sp.apodization.Gaussian(FWHM="100 Hz"),
            sp.FFT(),
        ]
    )
    hahn_dataset = processor.apply_operations(dataset=sim.methods[0].simulation)
    solid_dataset = processor.apply_operations(dataset=sim.methods[1].simulation)

.. skip: next

.. plot::
    :context: close-figs

    fig, ax = plt.subplots(1, 2, subplot_kw={"projection": "csdm"}, figsize=[8.5, 3])
    ax[0].set_title("Hahn-Echo Spectrum")
    ax[0].plot(hahn_dataset.real)
    ax[0].invert_xaxis()
    ax[0].grid()
    ax[1].set_title("Solid-Echo Spectrum")
    ax[1].plot(solid_dataset.real)
    ax[1].invert_xaxis()
    ax[1].grid()
    plt.tight_layout()
    plt.show()

In the Hahn-echo spectrum, the :math:`\text{p}_I`-dependent frequency
contributions (i.e., the shielding contributions) were averaged to zero, leaving
only the :math:`\text{d}_I`-dependent frequency contributions (i.e., the
first-order quadrupolar contribution). Conversely, in the solid-echo spectrum,
the :math:`\text{d}_I`-dependent frequency contributions (i.e., the first-order
quadrupolar contribution) were averaged to zero, leaving only the
:math:`\text{p}_I`-dependent frequency contributions (i.e., the shielding
contributions).

While these two examples nicely illustrate numerous important concepts for
building custom methods, it should also be noted that identical spectra could
have been obtained with a simpler custom method that used the ``freq_contrib``
to remove the undesired frequency contributions. The code for these two methods
is illustrated below.

.. plot::
    :context: close-figs

    quad_only = Method(
        channels=["2H"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            SpectralDimension(
                count=512,
                spectral_width=2e4,  # in Hz
                events=[
                    SpectralEvent(
                        transition_queries=[{"ch1": {"P": [-1]}}],
                        freq_contrib=["Quad1_2"]
                    )
                ],
            )
        ],
    )

    shielding_only = Method(
        channels=["2H"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            SpectralDimension(
                count=512,
                spectral_width=2e4,  # in Hz
                events=[
                    SpectralEvent(
                        transition_queries=[{"ch1": {"P": [-1]}}],
                        freq_contrib=["Shielding1_0", "Shielding1_2"],
                    )
                ],
            )
        ],
    )

    sim = Simulator()
    sim.spin_systems = [SpinSystem(sites=[deuterium])]
    sim.methods = [quad_only, shielding_only]
    sim.run()

    processor = sp.SignalProcessor(
        operations=[
            sp.IFFT(),
            sp.apodization.Gaussian(FWHM="100 Hz"),
            sp.FFT(),
        ]
    )
    quad_only_dataset = processor.apply_operations(dataset=sim.methods[0].simulation)
    shielding_only_dataset = processor.apply_operations(dataset=sim.methods[1].simulation)

.. skip: next

.. plot::
    :context: close-figs

    fig, ax = plt.subplots(1, 2, subplot_kw={"projection": "csdm"}, figsize=[8.5, 3])
    ax[0].set_title("Quad. Only Spectrum")
    ax[0].plot(quad_only_dataset.real)
    ax[0].invert_xaxis()
    ax[0].grid()
    ax[1].set_title("Shielding Only Spectrum")
    ax[1].plot(shielding_only_dataset.real)
    ax[1].invert_xaxis()
    ax[1].grid()
    plt.tight_layout()
    plt.show()



Origin and Reference Offset
---------------------------

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
to be no formal conventions for defining ``origin_offset`` and ``reference_offset``.


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
    - ``list``
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
    - ``list``
    - A list of :ref:`spectral_dim_api` objects describing the spectral dimensions for the method.

  * - affine_matrix
    - ``list`` or ``np.ndarray``
    - An *optional* (``n`` x ``n``) affine transformation matrix represented by a numpy array where ``n`` is
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
    - An *optional* float representing the origin offset, or Larmor frequency, along the
      spectroscopic dimension in units of Hz. The default value is ``None`` and the origin offset
      is set to the Larmor frequency of isotope from the :attr:`~mrsimulator.Method.channels`
      attribute of the method containing the spectral dimension.

  * - events
    - ``List``
    - An *optional* list of :ref:`event_api` used to emulate an experiment.
      The default value is an empty list.


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
      ``None`` and takes the global magnetic flux density defined by the method's
      :attr:`~mrsimulator.Method.magnetic_flux_density` attribute.

  * - rotor_angle
    - ``float``
    - An *optional* float describing the angle between the sample rotation axis and the external
      magnetic field in radians. The default is ``None`` and takes the global rotor angle defined
      by the method's :attr:`~mrsimulator.Method.rotor_angle` attribue.

  * - rotor_frequency
    - ``float``
    - An *optional* float describing the sample rotation frequency in Hz. For example, ``2000`` Hz.
      The default value is ``None`` and takes the global rotor frequency defined by the method's
      :attr:`~mrsimulator.Method.rotor_frequency` attribute.

  * - freq_contrib
    - ``List``
    - An *optional* list of :ref:`freq_contrib_api` (list of allowed strings) selecting which
      contributions to include when calculating a transition frequency. For example,
      ``["Shielding1_0", "Shielding1_2"]``. By default, the list is all frequency enumerations and
      all frequency contributions are calculated.

  * - transition_queries
    - ``list``
    - An *optional* ``list`` of :ref:`query_api` objects, or their ``dict`` representations
      selecting transitions active during the event. Only these selected transitions will
      contribute to the net frequency. The default is one TransitionQuery with *P=[0]*
      on ``ch1`` and ``None`` on all other channels


.. cssclass:: table-bordered table-striped centered
.. _table_mixing_event:
.. list-table:: The attributes of a MixingEvent object
  :widths: 20 15 65
  :header-rows: 1

  * - Attribute Name
    - Type
    - Description

  * - query
    - ``dict`` or :py:class:`~mrsimulator.method.MixingQuery`
    - A :py:class:`~mrsimulator.method.MixingQuery` object, or its ``dict`` representation,
      determines the complex amplitude of mixing between transitions in adjacent spectral
      or delay events.

..   - The coordinates along each spectral dimension are
..       described with the keywords,``count``(:math:`N`), ``spectral_width``
..       (:math:`\nu_\text{sw}`), and``reference_offset``(:math:`\nu_0`). The
..       coordinates are evaluated as,
..
..       .. math
..         \left([0, 1, 2, ... N-1] - \frac{T}{2}\right) \frac{\nu_\text{sw}}{N} + \nu_0
..
..       where :math:`T=N` when :math:`N` is even else :math:`T=N-1`.
