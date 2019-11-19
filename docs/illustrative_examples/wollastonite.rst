
Wollastonite
------------

Wollastonite is a high-temperature calcium-silicate,
:math:`\beta-\text{Ca}_3\text{Si}_3\text{O}_9`, with three distinct
:math:`^{29}\text{Si}` sites.  The :math:`^{29}\text{Si}` shielding tensor
parameters were obtained from
`Hansen et. al. <https://pubs.acs.org/doi/10.1021/ic020647f>`_


**Step 1**  Create sites.

.. doctest::

    >>> S29_1 = Site(isotope='29Si', isotropic_chemical_shift=-89.0, shielding_symmetric={'zeta': 59.8, 'eta': 0.62})
    >>> S29_2 = Site(isotope='29Si', isotropic_chemical_shift=-89.5, shielding_symmetric={'zeta': 52.1, 'eta': 0.68})
    >>> S29_3 = Site(isotope='29Si', isotropic_chemical_shift=-87.8, shielding_symmetric={'zeta': 69.4, 'eta': 0.60})

**Step 2**  Create isotopomers from these sites.

.. doctest::

    >>> isotopomers = [Isotopomer(sites=[site]) for site in [S29_1, S29_2, S29_3]]


**Step 3**  Create a dimension.

.. doctest::

    >>> dimension = Dimension(
    ...     isotope='29Si',
    ...     magnetic_flux_density=14.1, # in T
    ...     number_of_points=2046,
    ...     spectral_width=25000, # in Hz
    ...     reference_offset=-10000, # in Hz
    ...     rotor_frequency=1500 # in Hz
    ... )

**Step 4**  Create the Simulator object and add dimension and isotopomer
objects.

.. doctest::

    >>> sim_wollastonite = Simulator()
    >>> # add isotopomers
    >>> sim_wollastonite.isotopomers += isotopomers
    >>> # add dimensions
    >>> sim_wollastonite.dimensions += [dimension]

**Step 5**  Simulate and plot.

.. doctest::

    >>> x, y = sim_wollastonite.run(method=one_d_spectrum)
    >>> plot(x,y) # doctest:+SKIP

.. .. testsetup::
..     >>> plot_save(x, y, 'illustrative_example_wollastonite')

.. figure:: ../_images/illustrative_example_wollastonite.*
    :figclass: figure-polaroid
