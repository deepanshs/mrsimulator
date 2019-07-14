
.. _getting_started:

==================================
Getting started with `mrsimulator`
==================================

In this section, we have put together guidelines to encourage consistency
amongst the users. The NMR lineshapes are calculated with the instance of
the :ref:`simulator_api` class. Import this class following

.. doctest::

    >>> from mrsimulator import Simulator

and create an instance,

.. doctest::

    >>> sim1 = Simulator()

Here, ``sim1`` is an instance of the :ref:`simulator_api` class. The
two often used attributes of this instance are
attr:`~mrsimulator.Simulator.isotopomers` and
attr:`~mrsimulator.Simulator.spectrum`.
The default value of these attributes is

.. doctest::

    >>> print(sim1.isotopomers)
    []
    >>> print(sim1.spectrum)
    {}

Please read the section on :ref:`isotopomers` and :ref:`spectrum` for
instructions on how to assign values to these attributes, respectively.
Here, we use a built-in example to assign these attributes. Import
the examples using

.. doctest::

    >>> from mrsimulator import examples

and assign the isotopomers and spectrum attributes of ``sim1`` with

.. doctest::

    >>> sim1.isotopomers, sim1.spectrum = examples.csa_static()

You may view the contents of the isotopomers and spectrum attribute
using the print statement. Here, we use the pprint, the data pretty printer
command.

.. doctest::

    >>> from pprint import pprint
    >>> pprint(sim1.isotopomers)
    [{'abundance': '100 %',
      'sites': [{'isotope_symbol': '1H',
                 'isotropic_chemical_shift': '0 ppm',
                 'shielding_symmetric': {'anisotropy': '13.89 ppm',
                                         'asymmetry': 0.25}}]}]
    >>> pprint(sim1.spectrum)
    {'direct_dimension': {'magnetic_flux_density': '9.4 T',
                          'nucleus': '1H',
                          'number_of_points': 2048,
                          'reference_offset': '0 Hz',
                          'rotor_angle': '54.735 deg',
                          'rotor_frequency': '0 kHz',
                          'spectral_width': '25 kHz'}}

In general, the isotopomers contain the metadata on the spin system while
the spectrum contains metadata required to simulate the lineshapes.
A lineshape is simulated using the meth:`~mrsimulator.Simulator.run` method
of the :ref:`simulator_api` instance based on the NMR method.
In version 0.1, we provide `one_d_spectrum` method for simulating one
dimensional NMR lineshapes. Import this method using

.. doctest::

    >>> from mrsimulator.methods import one_d_spectrum

and run the simulation.

.. doctest::

    >>> freq, amp = sim1.run(one_d_spectrum, verbose=1)
    Setting up the virtual NMR spectrometer
    ---------------------------------------
    Adjusting the magnetic flux density to 9.4 T.
    Setting rotation angle to 0.9553059660790962 rad.
    Setting rotation frequency to 0.0 Hz.
    Detecting 1H(I=0.5, precession frequency = 400.228301848 MHz) isotope.
    Recording 1H spectrum with 2048 points over a 25000.0 Hz bandwidth and a reference offset of 0.0 Hz.
    <BLANKLINE>
    1H site 0 from isotopomer 0 @ 100.0% abundance
    ----------------------------------------------
    Isotropic chemical shift = 0.0 ppm
    Shielding anisotropy = 13.89 ppm
    Shielding asymmetry = 0.25

In the above code, the ``freq`` and ``amp`` are the frequency in Hz and the
corresponding amplitude of the spectrum. The following is a figure of the above
lineshape plotted using the matplotlib library.

.. doctest::

    >>> import matplotlib.pyplot as plt
    >>> def plot(x, y):
    ...     plt.plot(x,y)
    ...     plt.xlabel(f'frequency / {str(x.unit)}')
    ...     plt.show()

    >>> plot(freq, amp)

.. image:: /_static/1H_example.png
