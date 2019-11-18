
.. _config_simulator:

.. testsetup::

    >>> import matplotlib
    >>> font = {'family': 'Helvetica', 'weight': 'light', 'size': 9}
    >>> matplotlib.rc('font', **font)
    >>> from os import path
    >>> import matplotlib.pyplot as plt
    >>> def plot_save(x, y, filename):
    ...     plt.figure(figsize=(4, 3))
    ...     plt.plot(x, y, linewidth=1)
    ...     plt.xlim([x.value.max(), x.value.min()])
    ...     plt.xlabel(f"frequency ratio / {str(x.unit)}", **font)
    ...     plt.grid(color='gray', linestyle='--', linewidth=1.0, alpha=0.25)
    ...     plt.tight_layout(h_pad=0, w_pad=0, pad=0)
    ...
    ...     filename = path.split(filename)[1]
    ...     filepath = './docs/_images'
    ...     pth = path.join(filepath, filename)
    ...     plt.savefig(pth+'.pdf')
    ...     plt.savefig(pth+'.png', dpi=100)
    ...     plt.close()

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
    >>> sim.isotopomers = isotopomers
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
    The lineshape broadening seen in the above spectrum is the result of
    lineshapes arising from a distribution of tensors. In this case,
    the lineshape is an integral of 1000 individual spectrum.

Here is a spinning sideband spectrum at a 90-degree angle.

.. doctest::

    >>> sim.dimensions = [Dimension(isotope='29Si', spectral_width=25000, reference_offset=-7000, rotor_frequency=5000, rotor_angle=1.57079)]
    >>> x, y = sim.run(method=one_d_spectrum) # doctest:+SKIP
    >>> plt.plot(x,y) # doctest:+SKIP

.. .. testsetup::
..     >>> plot_save(x, y, 'example_amorphous_2')

.. figure:: _images/example_amorphous_2.*
    :figclass: figure-polaroid
