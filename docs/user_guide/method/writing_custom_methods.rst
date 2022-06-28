.. _writing_custom_methods:

======================
Writing Custom Methods
======================

The power of mrsimulator lies in the method object. Most time-domain NMR simulation
software requires users to walk through all aspects of a pulse sequence, but mrsimulator
preforms calculations in the frequency-domain and only needs descriptions of the spectral
dimensions being simulated. For a more in-depth discussion on how a method defined
spectral dimensions, see the :ref:`method_documentation`.

On this page, we will illustrate how to write custom methods by simulating progressively
more complicated NMR experiments on :math:`RbNO_3`. First lets import the necessary
classes and modules.

.. The Method object is where the versatility of mrsimulator becomes clear.
.. Most NMR density matrix simulations do all the calculations in the
.. time-domain, but mrsimulator performs its calculations in the frequency
.. domain. In these time-domain programs, you may set up an experiment that
.. walks through all aspects of a pulse sequence, but in mrsimulator, you
.. only need to set up a method describing all the spectral dimensions you
.. are simulating.
..
.. Each Method object holds global parameters, like magnetic_flux_density,
.. and a list of SpectralDimension objects, each one describing a dimension
.. of a multi-dimensional spectrum. Each SpectralDimension object contains
.. a list of events, in which you can adjust parameters, like rotor speed
.. or angle, select transitions based on their :math:`p` or :math:`d`
.. symmetries, etc. To illustrate this, let’s look at a few different
.. common NMR experiments on :math:`RbNO_3`, starting with a simple 1D
.. pulse-acquire experiment. We begin by making all necessary imports.

.. plot::
    :context: close-figs

    from mrsimulator import Site, SpinSystem, Simulator
    from mrsimulator.method import Method, SpectralDimension, SpectralEvent, MixingEvent
    from mrsimulator.spin_system.tensors import SymmetricTensor
    from mrsimulator import signal_processing as sp
    import matplotlib.pyplot as plt
    from pprint import pprint

Now we build the :math:`RbNO_3` spin system.

.. plot::
    :context: close-figs

    site1 = Site(
        isotope="87Rb",
        isotropic_chemical_shift=-27.4,  # ppm,
        quadrupolar=SymmetricTensor(Cq=1.68e6, eta=0.2)  # Cq in Hz
    )
    site2 = Site(
        isotope="87Rb",
        isotropic_chemical_shift=-28.5,  # ppm
        quadrupolar=SymmetricTensor(Cq=1.94e6, eta=1)  # Cq in Hz
    )
    site3 = Site(
        isotope="87Rb",
        isotropic_chemical_shift=-31.3,  # ppm
        quadrupolar=SymmetricTensor(Cq=1.72e6, eta=0.5)  # Cq in Hz
    )

    sites = [site1, site2, site3]
    spin_systems = [SpinSystem(sites=[s]) for s in sites]

One-Pulse Acquire
-----------------

Here we build our first method to simulate a simple one-pulse acquire experiment using the
generic Method object.

.. Now, we build the method. We will be building it from the generic Method
.. object, but you could just as easily use the built-in BlochDecaySpectrum
.. method.

.. plot::
    :context: close-figs

    pulseacquire = Method(
        channels=["87Rb"],
        magnetic_flux_density=9.4,  # in T
        rotor_frequency=10000,  # in Hz
        spectral_dimensions=[
            SpectralDimension(
                count=1e6,
                spectral_width=3e6,
                reference_offset=-3.5e3,
                events=[
                    SpectralEvent(transition_query=[{"ch1": {"P": [-1]}}])
                ]
            )

        ]
    )

The *channels* key holds the nucleus being probed, here rubidium-87. The
*magnetic_flux_density* key holds the external magnetic field stength in T, and
*rotor_frequency* holds the rotor frequency in Hz. The *spectral_dimensions* key
holds a list of SpectralDimension objects defining the spectral grid of the method.
each SpectralDimension object contains a *count* key, defining the number of points
in that dimension, a *spectral_width* key, containing the spectral width in Hz,
and a *reference_offset* key, containing the reference offset of the dimension in Hz.
Each spectral dimension also contains an *events* key which is a list events defining
the coherences to simulate. In this example, we are using a spectral event to
select the :math:`p=m_f-m_i=-1` coherences for this simulation.

Next we set up and run the simulator object with our spin system and method.

.. plot::
    :context: close-figs
    :caption: A simulated one-pulse acquire spectrum of :math:`{87}^\text{Rb}` with all sidebands shown (left) and zoomed in plot of the central transition (right).

    sim = Simulator()
    sim.spin_systems = spin_systems
    sim.methods = [pulseacquire]
    sim.config.number_of_sidebands = 256
    sim.run()

Now, we create a signal processing object to add some exponential line broadening
to the simulated spectrum and plot the processed dataset.

.. skip: next

.. plot::
    :context: close-figs

    processor = sp.SignalProcessor(
        operations=[
            sp.IFFT(),
            sp.apodization.Exponential(FWHM="10 Hz"),
            sp.FFT(),
        ]
    )

    processed_data = processor.apply_operations(data=sim.methods[0].simulation.real)

    fig, ax = plt.subplots(
        nrows=1,
        ncols=2,
        subplot_kw={"projection": "csdm"},
        figsize=(8, 4)
    )

    ax[0].plot(processed_data.real, color="black", linewidth=1)
    ax[0].invert_xaxis()
    ax[1].plot(processed_data.real, color="black", linewidth=1)
    ax[1].set_xlim(-50, 0)
    ax[1].invert_xaxis()
    plt.tight_layout()
    plt.show()

Selecting the Central Transition
--------------------------------

Now, let’s say we wanted to supress the satellites. To do this, we need
to simulate a central-transition-selective 1D experiment. We now add a restriction to
:math:`D`, defined as :math:`D = m_f^2 -m_i^2`, in our transition query. For the
central-transition selective method, we specify :math:`D=0`.

.. plot::
    :context: close-figs

    ct_pulseacquire = Method(
        channels=["87Rb"],
        magnetic_flux_density=9.4,  # in T
        rotor_frequency=10000,  # in Hz
        spectral_dimensions=[
            SpectralDimension(
                count=20000,
                spectral_width=8e3,
                reference_offset=-3.5e3,
                events=[
                    SpectralEvent(transition_query=[{"ch1": {"P": [-1], "D": [0]}}])
                ]
            )
        ]
    )

We now replace the old ``pulseacquire`` method in the simulator object with our new
``ct_pulseacquire`` method and re-simulate the spectrum.

.. We simply add this new method to the simulator object, run the
.. simulation, apply our proceessing, and plot the data.

.. skip: next

.. plot::
    :context: close-figs
    :caption: A simulated central-transition selective spectrum of :math:`{87}^\text{Rb}`. The large number of sidebands from the previous simulation have been suppressed.

    sim.methods = [ct_pulseacquire]
    sim.config.number_of_sidebands = 70  # Reset number of sidebands for efficiency
    sim.run()

    processed_data = processor.apply_operations(data=sim.methods[0].simulation.real)

    plt.figure(figsize=(6, 4))
    ax = plt.subplot(projection="csdm")
    # ax.plot(sim.methods[0].simulation, color="blue", linewidth=1)
    ax.plot(processed_data.real, color="black", linewidth=1)
    ax.invert_xaxis()
    plt.tight_layout()
    plt.show()

Three-Quantum MAS
-----------------

Now, let’s construct a method to simulate a 3Q-MAS spectrum.

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
                # label="Isotropic dimension",
                events=[
                    SpectralEvent(transition_query=[{"ch1": {"P": [-3], "D": [0]}}])
                ]
            ),
            SpectralDimension(
                count=256,
                spectral_width=6e3,  # in Hz
                reference_offset=-5e3,  # in Hz
                # label="MAS dimension",
                events=[
                    SpectralEvent(transition_query=[{"ch1": {"P":[-1], "D": [0]}}])
                ]
            )
        ],
    )

Now, instead of just one item in the list of spectral dimensions, we
have two, because 3Q-MAS is a two-dimensional experiment. In the first
dimension, we are selecting the triple-quantum coherence by specifying a
transition query of :math:`p=-3` and :math:`d=0`. In the MAS dimension,
we are selecting the central transition with a transition query of
:math:`p=-1` and :math:`d=0`.

Again, we add this method to the simulator object, run the simulation, and
plot the data.

.. skip: next

.. plot::
    :context: close-figs
    :caption: An unsheared 3Q-MAS spectrum of :math:`{87}^\text{Rb}`

    sim.methods = [mqmas]
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
    data = processor.apply_operations(data=sim.methods[0].simulation)

    plt.figure(figsize=(6, 4))
    ax = plt.subplot(projection="csdm")
    cb = ax.imshow(data.real / data.real.max(), aspect="auto", cmap="gist_ncar_r")
    plt.colorbar(cb)
    ax.invert_xaxis()
    ax.invert_yaxis()
    plt.tight_layout()
    plt.show()

Sheared Three-Quantum MAS
-------------------------

For 3Q-MAS experiments, however, the spectrum is often sheared and
scaled to make the vertical dimension the purely isotropic dimension.
This can be accomplished with an affine matrix added to the method.
Let’s re-make our 3Q-MAS method with this affine matrix.

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

.. note:
    The *affine_matrix* in mrsimulator is given in row-major as a n by n array
    where n is the number of spectral dimensions

Again, we now add the method to the simulator object, run the
simulation, and plot the data.

.. skip: next

.. plot::
    :context: close-figs
    :caption: A 3Q-MAS spectrum of :math:`{87}^\text{Rb}` sheared such that the dimensions are purely MAS and isotropic.

    sim.methods = [sheared_mqmas]
    sim.run()

    data = processor.apply_operations(data=sim.methods[0].simulation)

    plt.figure(figsize=(6, 4))
    ax = plt.subplot(projection="csdm")
    cb = ax.imshow(data.real / data.real.max(), aspect="auto", cmap="gist_ncar_r")
    plt.colorbar(cb)
    ax.set_ylim((-70, -47))
    ax.invert_xaxis()
    ax.invert_yaxis()
    plt.tight_layout()
    plt.show()


For convenience sake, the one-pulse acquire (BlochDecaySpectrum), one-pulse acquire central
transition selective (BlochDecayCTSpectrum), and Three-Quantum MAS (ThreeQ_VAS) methods
along with other common methods can be imported from the ``mrsimulator.method.lib`` package.
For more details, see the :ref:`methods_library_documentation`.

.. For the convenience methods mentioned here and more, please see our
.. methods library. For a more in-depth description of creating methods,
.. see our advanced users methods page.

Hahn vs Solid Echo
------------------

We have seen how a Method object can select between different coherences by using
SpectralDimension and SpectralEvents. By adding a MixingEvent, we can selectively simulate
frequencies from specific transition pathways. Below we construct a deuterium spin system
and two Method objects to simulate a Hahn and Solid Echo experiment.

.. plot::
    :context: close-figs

    deuterium = Site(
        isotope="2H",
        isotropic_chemical_shift=10,  # in ppm
        shielding_symmetric=SymmetricTensor(zeta=-80, eta=0.25),  # zeta in ppm
        quadrupolar=SymmetricTensor(Cq=10e3, eta=0.0))

    spin_system = SpinSystem(sites=[deuterium])

Hahn Echo
"""""""""

The Hahn Echo experiment observes the transition frequencies from the following
:math:`\mathbb{p}` transition symmetry pathways (a.k.a coherence transfer pathways).
For more discussion on transition symmetry pathways, see the ((pathway documentation page??))

.. math::

    \mathbb{p}: 0 \xrightarrow[]{\frac{\pi}{2}} +1 \xrightarrow[]{\pi} -1

(??) This pathway selectively refocuses the :math:`\mathbb{p}` frequency contributions into
an echo while leaving the :math:`\mathbb{d}` contributions free to evolve unaffected by the
:math:`\pi` pulse. (??)
Below is a diagram representing the different energy level transitions and corresponding
pathways observed by the Hahn Echo experiment.

.. figure:: ../../_static/deuteriumHahnEcho.*
    :alt: Transition symmetry pathways for the Hahn Echo experiment
    :align: center
    :width: 50%

    Energy level transitions and symmetry pathways for the Hahn Echo experiment.

Although a normal experiment would start with a :math:`\frac{\pi}{2}` rotation to transfer the
equilibrium magnetization to a desired symmetry, we can eliminate the first rotation in
mrsimulator by defining the first symmetry as :math:`\mathbb{p} = +1`. Our transition symmetry
pathway now becomes

.. math::

    \mathbb{p}: +1 \xrightarrow[]{\pi} -1

Below is a method object which simulated the Hahn Echo experiment. The MixingEvent defines the
:math:`\pi` rotation between the two SpectralEvents. We also plot the transition pathways for

.. plot::
    :context: close-figs

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

    pprint(hahn_echo.get_transition_pathways(spin_system))

.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    [|1.0⟩⟨0.0| ⟶ |-1.0⟩⟨0.0|, weight=(1+0j)
     |0.0⟩⟨-1.0| ⟶ |0.0⟩⟨1.0|, weight=(1+0j)]



Solid Echo
""""""""""

Any discussion such as transition pathways goes here

.. figure:: ../../_static/deuteriumSolidEcho.*
    :alt: Transition symmetry pathways for the Hahn Echo experiment
    :align: center
    :width: 50%

    Energy level transitions and symmetry pathways for the Solid Echo experiment.

.. math::

    \mathbb{p}: 0 \xrightarrow[]{\frac{\pi}{2}} +1 \xrightarrow[]{\frac{\pi}{2}} -1

simplifies to

.. math::

    \mathbb{p}: -1 \xrightarrow[]{\frac{\pi}{2}} -1

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



Now we setup and run the simulation then process and plot the data

.. skip: next

.. plot::
    :context: close-figs
    :caption: Simulated Hanh Echo spectrum (left) and Solid Echo spectrum (right) for the same :math:`2^\text{H}` spin system.

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
    hahn_data = processor.apply_operations(data=sim.methods[0].simulation)
    solid_data = processor.apply_operations(data=sim.methods[1].simulation)

    fig, ax = plt.subplots(
        nrows=1,
        ncols=2,
        subplot_kw={"projection": "csdm"},
        figsize=[8, 4]
    )

    ax[0].plot(hahn_data.real, color="black", linewidth=1)
    ax[0].invert_xaxis()
    ax[1].plot(solid_data.real, color="black", linewidth=1)
    ax[1].invert_xaxis()
    plt.tight_layout()
    plt.show()
