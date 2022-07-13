.. _transition_query_documentation:

==================
Transition Queries
==================

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
**SpectralEvent** or **DelayEvent**. Keep in mind, however, that **Method**
objects are designed without any details of the spin systems upon which they
will act. For example, in the density matrix of a spin system ensemble, one
could easily identify a transition by its row and column indexes. However, those
indexes depend on the spin system and how the spins and their eigenstates have
been assigned to those indexes. Instead, we need a spin-system agnostic approach
for selecting transitions.


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

.. figure:: ../_static/SpinOneThreeHalves.*
    :width: 500
    :alt: figure
    :align: center

.. figure:: ../_static/SpinFiveHalf.*
    :width: 550
    :alt: figure
    :align: center

----

For discussion on higher-level transition symmetry functions and how they can be used when
selecting NMR coherence order pathways, click on the following dropdown.

.. raw:: html

    <details>
    <summary><b>Higher-Level Symmetry Functions</b></summary>

.. include:: spin_trans_symm.rst

.. raw:: html

    </details>

Single-Spin Queries
-------------------

Based on the review above, we now know for the spin :math:`I=1`, the transition
:math:`\ketbra{-1}{0}` can be selected with :math:`{(\text{p}_I,\text{d}_I) =
(-1,1)}`.  In **mrsimulator**, this transition is selected during a
**SpectralEvent** using the **SymmetryQuery** and **TransitionQuery** objects,
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
    passed to a higher level object. Since a **SpectralEvent** object holds a list of
    **TransitionQuery** objects, the above code could have been written as

    .. plot::
        :context: close-figs

            # Both SymmetryQuery and TransitionQuery object as dict
            symm_query_dict = {"P": [-1], "D": [1]}
            trans_query_dict = {"ch1": symm_query_dict}

            # Dictionary of TransitionQuery passed to SpectralEvent
            spec_event = SpectralEvent(transition_queries=[trans_query_dict])

In the example above, the **SymmetryQuery** object is created and assigned to
the **TransitionQuery** attribute ``ch1``, i.e., it acts on the isotope in the
"first channel". Recall that the ``channels`` attribute of the **Method** object
holds an ordered list of isotope strings. This list's first, second, and third
isotopes are associated with ``ch1``, ``ch2``, and ``ch3``, respectively.
Currently, **mrsimulator** only supports up to three channels, although this may
be increased in future versions.

The **TransitionQuery** object goes into an unordered list in the
``transition_queries`` attribute of a **SpectralEvent** object. The
**SpectralEvent** object, in turn, is added to an ordered list in the ``events``
attribute of a **SpectralDimension** object. All this is illustrated in the code
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

    # This method selects all observable (p_I=–1) transitions
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

    # This method selects observable (p_I=–1) transitions with d_I = 1
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

    # This method selects observable (p_I=–1) transitions with d_I = -1
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

    Whenever the ``D`` attribute is omitted, the **SymmetryQuery** allows
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
**SpectralDimension** objects, each holding a single **SpectralEvent**.  The
**TransitionQuery** objects select each transition in their respective
**SpectralDimension** objects.

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
to map the dataset amplitude to colors; The `colorbar()
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


.. **Method** also has a related function :py:meth:`~mrsimulator.Method.plot` for generating a symmetry
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

.. figure:: ../_static/ThreeCoupledSpinsEnergy.*
    :width: 550
    :alt: figure
    :align: center

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
    ax.plot(processor.apply_operations(dataset=sim.methods[0].simulation))
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
"indistinguishable" to a **TransitionQuery** object. Recall that ``ch1`` is
associated with the first isotope in the list of isotope strings assigned to the
**Method** attribute ``channels``. When the **TransitionQuery** above is combined with
the **SpinSystem** object with three :math:`^1\text{H}` Sites, it must first expand
its **SymmetryQuery** into an intermediate set of spin-system-specific symmetry
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
you to inspect the transitions selected by the **Method**
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

.. figure:: ../_static/ThreeCoupledSpinsDoubleQuantum.*
    :width: 500
    :alt: figure
    :align: center

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

.. figure:: ../_static/ThreeCoupledSpinsSingleQuantum.*
    :width: 500
    :alt: figure
    :align: center

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
    ax.plot(processor.apply_operations(dataset=sim.methods[0].simulation))
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

.. figure:: ../_static/Spin1SpinHalfCouple.*
    :width: 850
    :alt: figure
    :align: center

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
be obtained through the use of multiple **SpectralEvent** objects in  the
**SpectralDimension** associated with the isotropic dimension, as shown in the
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
``affine_matrix`` in the **Method** object to do this shear is given by

``affine_matrix=[[1,0],[-8/25, 17/25]]``
