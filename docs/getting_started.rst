
.. _getting_started:

==================================
Getting started with `mrsimulator`
==================================

The NMR lineshape is calculated through the instance of the
:ref:`simulator_api` class. Import the class

.. doctest::

    >>> from mrsimulator import Simulator

and create an instance,

.. doctest::

    >>> sim1 = Simulator()

Here, ``sim1`` is an instance of the :ref:`simulator_api` class. The
two often used attributes of this instance are
:attr:`~mrsimulator.Simulator.isotopomers` and
:attr:`~mrsimulator.Simulator.spectrum`.
The default value of these attributes is

.. doctest::

    >>> print(sim1.isotopomers)
    []
    >>> print(sim1.spectrum)
    {}

Read the section on :ref:`isotopomers` and :ref:`spectrum` for instructions
on how to assign the values to these attributes, respectively.
Here, we use a built-in example to assign these attributes.

.. doctest::

    >>> from pprint import pprint
    >>> from mrsimulator import examples
    >>> sim1.isotopomers, sim1.spectrum = examples.csa_static()
    >>> pprint(sim1.isotopomers)
    [{'abundance': '100 %',
      'sites': [{'isotope_symbol': '1H',
                 'isotropic_chemical_shift': '0 Hz',
                 'shielding_symmetric': {'anisotropy': '-3.89 kHz',
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
To simulate a lineshape, use the :meth:`~mrsimulator.Simulator.run` method
of the instance with an NMR method.
In version 0.1, there is only one method, `one_d_spectrum`.
Import this method using

.. doctest::

    >>> from mrsimulator.methods import one_d_spectrum

and run the simulation.

.. doctest::
    :skipif: None is None

    >>> freq, amp = sim1.run(one_d_spectrum, verbose=1)
    <BLANKLINE>
    Setting up the virtual NMR spectrometer
    ---------------------------------------
    Adjusting the magnetic flux density to 9.4 T
    Setting rotation angle to 0.9553059660790962 rad
    Setting rotation frequency to 0.0 Hz
    Detecting 1H(I=0.5, precession frequency = 400.228301848 MHz) isotope
    Recording 1H spectrum with 2048 points over a 25000.0 Hz bandwidth and a reference offset of 0.0 Hz.
    <BLANKLINE>
    1H site 0 in isotopomer 0
    ----------------------------
    isotropic chemical shift = 0.0 Hz
    chemical shift anisotropy = -3890.0 Hz
    chemical shift asymmetry = 0.25
    <BLANKLINE>
    Execution time 0.015989 s

In the above code, the ``freq`` and ``amp`` are the frequency and the
corresponding amplitude of the spectrum. The following is a plot of the
lineshape using the matplotlib library.

..doctest::

    >>> import matplotlib.pyplot as plt
    >>> plt.plot(freq, amp)
    >>> plt.xlabel('frequency / Hz')
    >>> plt.show()

.. image:: /_static/1H_example.pdf
