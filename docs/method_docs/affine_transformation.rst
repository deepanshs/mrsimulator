.. _affine_transformation_documentation:

======================
Affine Transformations
======================

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

    import numpy as np
    import mrsimulator.signal_processor as sp
    from mrsimulator import Site, SpinSystem, Simulator
    from mrsimulator.method import Method, SpectralDimension, SpectralEvent

    # Create same three sites in RbNO3 from previous section
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

    gauss_convolve = sp.SignalProcessor(
        operations=[
            sp.IFFT(dim_index=(0, 1)),
            sp.apodization.Gaussian(FWHM="0.08 kHz", dim_index=0),
            sp.apodization.Gaussian(FWHM="0.22 kHz", dim_index=1),
            sp.FFT(dim_index=(0, 1)),
        ]
    )
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
