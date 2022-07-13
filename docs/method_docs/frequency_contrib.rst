.. _frequency_contrib_documentation:

=======================
Frequency Contributions
=======================

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
is the root of the second-rank Legendre polynomial ":math:`P_2(\cos \theta_R)`.
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

As described in ":ref:`method_theory`", these transition symmetry functions play an
essential role in evaluating the individual frequency contributions to the
overall transition frequency, given in the table below and in
:py:meth:`~mrsimulator.method.frequency_contrib.FrequencyEnum`. They also aid in
pulse sequence design by identifying how different frequency contributions
refocus through the transition pathways.

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

.. note::

    **Echo Symmetry Classification**

    The well-known Hahn-echo can occur whenever the :math:`p_I` values of
    transitions in a transition pathway change sign.  This is because the
    changing sign of :math:`p_I` leads to a sign change for every
    :math:`p_I`-dependent transition frequency contribution. Thus, a Hahn
    echo forms whenever

    .. math::
        \overline{\text{p}_I} = \frac{1}{t} \int_0^t \text{p}_I(t') \, dt' = 0,

    assuming a frequency contribution's spatial symmetry function, :math:`{\Xi}`,
    remains constant during this period.  As seen in the table in the
    :ref:`frequency_contribution_table` table, sign changes in other symmetry
    functions can also lead to corresponding sign changes for dependent
    frequency contributions.  Thus, a problem with showing only the :math:`p_I`
    symmetry pathway for an NMR method is that it does not explain the formation
    of other classes of echoes that result when other symmetry functions change
    sign in a transition pathway.  To fully understand when and which frequency
    contributions refocus into echoes, we must follow *all* relevant spatial,
    transition, or spatial-transition product symmetries through an NMR
    experiment.   Thus, we generally classify echoes that refocus during a time
    interval as a *transition symmetry echo* (at constant :math:`{\Xi}_k`) when

    .. math::
        \overline{{\xi}_k} = \frac{1}{t} \int_0^t {\xi}_k(t') \, dt' = 0,

    and as a *spatial symmetry echo* (at constant :math:`{\xi}_k`) when

    .. math::
        \overline{{\Xi}_k} = \frac{1}{t} \int_0^t {\Xi}_k(t') \, dt' = 0,

    and as a *spatial-transition symmetry product* echo when

    .. math::

        \overline{{\Xi}_k {\xi}_k} = \frac{1}{t} \int_0^t {\Xi}_k(t') \, {\xi}_k(t')  \, dt' = 0.

    Within the class of transition echoes we find subclasses such as
    :math:`\text{p}` echoes, which include the Hahn echo and the stimulated
    echo; :math:`\text{d}` echoes, which include the solid echo and Solomon
    echoes,  :math:`\text{c}_4` echoes, used in MQ-MAS and Satellite-Transition
    Magic-Angle Spinning (ST-MAS); :math:`\text{c}_2` echoes, used in
    Correlation Of Anisotropies Separated Through Echo Refocusing (COASTER); and
    :math:`\text{c}_0` echoes, used in Multiple-Quantum DOuble Rotation
    (MQ-DOR).

    Within the class of spatial echoes we find subclasses such as :math:`\mathbb{D}`
    rotary echoes, which occur during sample rotation, and :math:`\mathbb{D}_0` and
    :math:`\mathbb{G}_0` echoes, which are designed to occur simultaneously during the
    Dynamic-Angle Spinning (DAS) experiment.

p and d Echoes on Deuterium
'''''''''''''''''''''''''''

Here, we examine two examples in a deuterium spin system that illustrate the
importance of echo classification in understanding how transition-frequency
contributions can be eliminated or separated based on their dependence on
different transition symmetry functions.

First, we implement two **Method** objects that follow the design of the
experimental pulse sequence. In this effort, we use **RotationQuery** objects to
select the desired transition pathways and obtain spectra with the desired
average frequencies. Then, we implement two simpler **Method** objects that
produce identical spectra and illustrate how :ref:`frequency
contributions<freq_contrib_api>` can be used to reduce the number of events
needed in a custom method.

Consider the Hahn and Solid Echo pulse sequences on the left and right,
respectively.

.. figure:: ../_static/HahnAndSolidEcho.*
    :alt: Transition symmetry pathways for the Hahn and Solid Echo experiments
    :align: center
    :width: 100%

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

Below are two custom **Method** objects for simulating the Hahn and Solid Echo
experiments. There is only one **SpectralDimension** object in each method, and
the average frequency during each spectral dimension is derived from equal
fractions of two **SpectralEvent** objects.  Between these two **SpectralEvent**
objects is a **MixingEvent** with a **RotationQuery** object. The
**RotationQuery** object is created with a :math:`\pi` rotation in the Hahn Echo
method, and a :math:`\pi/2` rotation in the Solid Echo method.

.. note ::

    The ``transition_queries`` attribute of **SpectralEvent** holds a list of
    **TransitionQuery** objects. Each **TransitionQuery** in the list applies to
    the full set of transitions in the spin system. The union of these transition
    subsets becomes the final set of selected transitions during the
    **SpectralEvent**.

We use the deuterium Site defined earlier in this document.

.. plot::
    :context: close-figs

    import numpy as np
    from mrsimulator import Site, SpinSystem, Simulator
    from mrsimulator.method import Method, SpectralDimension, SpectralEvent, MixingEvent

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

    import mrsimulator.signal_processor as sp

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
