
Potassium Sulfate
-----------------

The following example is a :math:`^{33}\text{S}` NMR lineshape simulation of
potassium sulfate (:math:`\text{K}_2\text{SO}_4`). The quadrupole tensor
parameters for :math:`^{33}\text{S}` is obtained from
`Moudrakovski et. al. <https://pubs.acs.org/doi/10.1021/jp908206c>`_


**Step 1**  Create sites, in this case, just the one.

.. doctest::

    >>> S33 = Site(
    ...     name='33S',
    ...     isotope='33S',
    ...     isotropic_chemical_shift=335.7,
    ...     quadrupolar={'Cq': 0.959e6, 'eta': 0.42}
    ... )

**Step 2**  Create isotopomers from this site.

.. doctest::

    >>> isotopomers = Isotopomer(sites=[S33])


**Step 3**  Create a dimension.

.. doctest::

    >>> dimension = Dimension(
    ...     isotope='33S',
    ...     magnetic_flux_density=21.14, # in T
    ...     number_of_points=2046,
    ...     spectral_width=5000, # in Hz
    ...     reference_offset=22500, # in Hz
    ...     rotor_frequency=14000 # in Hz
    ... )

**Step 4**  Create the Simulator object and add dimension and isotopomer
object.

.. doctest::

    >>> sim_K2SO3 = Simulator()
    >>> sim_K2SO3.isotopomers += [isotopomers]
    >>> sim_K2SO3.dimensions += [dimension]

**Step 5**  Simulate and plot.

.. doctest::

    >>> x, y = sim_K2SO3.run(method=one_d_spectrum)
    >>> plt.plot(x,y) # doctest:+SKIP

.. .. testsetup::
..     >>> plot_save(x, y, 'illustrative_example_K2SO3')

.. figure:: ../_images/illustrative_example_K2SO3.*
    :figclass: figure-polaroid
