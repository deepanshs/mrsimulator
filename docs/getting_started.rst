

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

    >>> print(sim0.isotopomers)
    []
    >>> print(sim0.spectrum)
    {}

Read the section on :ref:`isotopomers` and :ref:`spectrum` for instructions
on how to assign the values to these attributes, respectively.
Here, we use a built-in example to assign these attributes.

.. doctest::

    >>> from mrsimulator import examples
    >>> sim1.isotopomers, sim1.spectrum = examples.csa_static()
    >>> print(sim1.isotopomers)
    [
      {
        "sites": [
          {
            "isotope_symbol": "1H",
            "isotropic_chemical_shift": "0 Hz",
            "shielding_symmetric": {
              "anisotropy": "-3.89 kHz",
              "asymmetry": 0.25
            }
          }
        ],
        "abundance": "100 %"
      }
    ]
    >>> print(sim1.spectrum)
    {
      "direct_dimension": {
        "magnetic_flux_density": "9.4 T",
        "rotor_frequency": "0 kHz",
        "rotor_angle": "54.735 deg",
        "number_of_points": 2048,
        "spectral_width": "25 kHz",
        "reference_offset": "0 Hz",
        "nucleus": "1H"
      }
    }

In general, the isotopomers contain the metadata on the spin system while
the spectrum contains metadata required to simulate the lineshapes.
To simulate a lineshape, use the :meth:`~mrsimulator.Simulator.run` method
of the instance with an NMR method.
In version 0.1, there is only one method, `one_d_spectrum`.
Import this method using

.. doctest::

    >>> from mrsimulator.methods import one_D_spectrum

and run the simulation.

.. doctest::

    >>> freq, amp = sim0.run(one_d_spectrum, verbose=1)
    Setting up the virtual NMR spectrometer
    ---------------------------------------
    Adjusting the magnetic flux density to 9.4 T
    Setting rotation angle to 0.9553059660790962 rad
    Setting rotation frequency to 0.0 Hz
    Detecting 1H(I=0.5, precession frequency = 400.228301848 MHz) isotope
    Recording 1H spectrum with 2048 points over a 25000.0 Hz bandwidth and a reference offset of 0.0 Hz.

    1H site 0 in isotopomer 0
    ----------------------------
    isotropic chemical shift = 0.0 Hz
    chemical shift anisotropy = -3890.0 Hz
    chemical shift asymmetry = 0.25

    Execution time 0.015989 s

In the above code, the ``freq`` and ``amp`` are the frequency and the
corresponding amplitude of the spectrum. The following is a plot of the
lineshape using the matplotlib library.

    >>> import matplotlib.pyplot as plt
    >>> plt.plot(freq, amp)
    >>> plt.xlabel('frequency / Hz')
    >>> plt.show()

.. image:: /_static/1H_example.pdf

