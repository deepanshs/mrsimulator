
.. _getting_started:

==================================
Getting started with `mrsimulator`
==================================

We have put together a set of guidelines from using various methods and
attributes of `mrsimulator` package. We encourage the users
to follow these guidelines to promote consistency amongst others.
The solid state nuclear magnetic resonance (ssNMR) lineshapes are calculated
through an instance of the :ref:`simulator_api` class.

First, import the :ref:`simulator_api` class,

.. doctest::

    >>> from mrsimulator import Simulator

and create an instance as follows,

.. doctest::

    >>> sim1 = Simulator()

Here, ``sim1`` is a variable with an instance of the :ref:`simulator_api`
class. The two attributes of this class that you will often use are
:attr:`~mrsimulator.Simulator.isotopomers` and
:attr:`~mrsimulator.Simulator.spectrum`. The default value of these
attributes is an empty list.

.. doctest::

    >>> sim1.isotopomers
    []
    >>> sim1.spectrum
    []

Before we can start simulating lineshapes, we need to specify a list of
isotopomers and spectrum to the above attributes.

Setting up an Isotopomer object
-------------------------------
We define an isotopomer as an isolated spin-system containing
multiple sites and couplings between them. In this example, however, we
concern ourselves with a single site spin-system, that is, an isotopomer with
a single site. Shown below is an example of such isotopomer expressed as a
python dictionary.

.. code-block:: json

    {
        "sites": [
            {
                "isotope": "29Si",
                "isotropic_chemical_shift": "-101.1 ppm",
                "shielding_symmetric": {
                    "anisotropy": "70.5 ppm",
                    "asymmetry": 0.1
            }
        ]
    }

The python dictionary contains a ``sites`` keyword containing a list of sites
defined within the isotopomer. In this examples, we have defined a single site,
again as a python dictionary, containing site specific information such as,
the site isotope, the isotropic chemical shift, and the parameters from the
irreducible second rank symmetric nuclear shielding tensor---asymmetry and
anisotropy, expressed in Haeberlen convension. For additional information
see :ref:`isotopomer` and :ref:`site`.




An instance of an Isotopomer object may be created from a python dictionary,
such as the one shown above, using the :ref:`isotopomer_api`
class.

    >>> from mrsimulator import Isotopomer
    >>> isotopomer_object = Isotopomer.parse_json_with_units(isotopomer_dict)

You may create as many isotopomers as necessary, although in this example, we
stick with a single isotopomer. Now add this isotopomer to the variable
``sim``



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
      'sites': [{'isotope': '1H',
                 'isotropic_chemical_shift': '0 ppm',
                 'shielding_symmetric': {'anisotropy': '13.89 ppm',
                                         'asymmetry': 0.25}}]}]
    >>> pprint(sim1.spectrum)
    {'direct_dimension': {'magnetic_flux_density': '9.4 T',
                          'isotope': '1H',
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
