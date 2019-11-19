
.. _amorphous_materials:

===================================
Simulating amorphous-like materials
===================================

In this section, we illustrate how ``Mrsimulator`` may be used in simulating
amorphous materials. We do this by assuming a distribution of tensors. For
example, consider the following,

.. doctest::

    >>> import numpy as np
    >>> iso = np.random.normal(loc=-100.0, scale=5.0, size=1000)
    >>> zeta = np.random.normal(loc=50.0, scale=12.12, size=1000)
    >>> eta = np.random.normal(loc=0.5, scale=0.1, size=1000)

Here, we have created three Gaussian distributions, one for each isotropic
chemical shift, ``iso``, shielding anisotropy, ``zeta``, and shielding
asymmetry, ``eta``. The isotropic chemical shift distribution is centered at
-100 ppm with a standard deviation of 5 ppm. Similarly, the ``zeta`` and
``eta`` distributions are centered at 50 ppm and 0.5 with 12.12 ppm and 0.1
standard deviations, respectively. A total of 1000 points is sampled from the
three distributions.

Let's create the site and isotopomers objects from these parameters.

.. doctest::

    >>> from mrsimulator import Simulator, Site, Isotopomer, Dimension
    >>> from mrsimulator.methods import one_d_spectrum
    >>> isotopomers = []
    >>> for i, z, e in zip(iso, zeta, eta):
    ...     site = Site(isotope='29Si', isotropic_chemical_shift=i, shielding_symmetric={'zeta':z, 'eta':e})
    ...     isotopomers.append(Isotopomer(sites=[site]))


Now, that we have a 1000 isotopomers, let's create the Simulator object and add
these isotopomers.

.. doctest::

    >>> sim = Simulator()
    >>> # add isotopomers
    >>> sim.isotopomers = isotopomers
    >>> # create and add a dimension
    >>> sim.dimensions = [Dimension(isotope='29Si', spectral_width=25000, reference_offset=-7000)]

Let's observe the static spectrum first.

.. doctest::

    >>> x, y = sim.run(method=one_d_spectrum) # doctest:+SKIP
    >>> plt.plot(x,y) # doctest:+SKIP

.. .. testsetup::
..     >>> plot_save(x, y, 'example_amorphous_1')

.. figure:: _images/example_amorphous_1.*
    :figclass: figure-polaroid

.. note::
    The broad lineshape seen in the above spectrum is the result of
    lineshapes arising from a distribution of tensors. In this case,
    the lineshape is an integral of 1000 individual spectrum. There is no
    lineshape broadening filter applied to the spectrum.

Here is another example with a spinning sideband spectrum at a 90-degree angle.

.. doctest::

    >>> sim.dimensions[0].rotor_frequency = 5000 # in Hz
    >>> sim.dimensions[0].rotor_angle = 1.57079 # 90 degree in radian
    >>> x, y = sim.run(method=one_d_spectrum) # doctest:+SKIP
    >>> plt.plot(x,y) # doctest:+SKIP

.. .. testsetup::
..     >>> plot_save(x, y, 'example_amorphous_2')

.. figure:: _images/example_amorphous_2.*
    :figclass: figure-polaroid
