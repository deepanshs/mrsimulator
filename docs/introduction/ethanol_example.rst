.. _introduction_ethanol_example:

Ethanol Example
^^^^^^^^^^^^^^^

Here we work through an example that should be familiar to nearly all practitioners of NMR spectroscopy, i.e., 
the simulation of the :math:`^1\text{H}` and :math:`^{13}\text{C}` liquid-state NMR spectra 
of ethanol with its various isotopomers. The :math:`^1\text{H}` spectrum will include the characteristic
:math:`^{13}\text{C}` `satellite peaks <https://en.wikipedia.org/wiki/Carbon-13_NMR_satellite>`_
which arise from couplings between :math:`^{1}\text{H}` and :math:`^{13}\text{C}` in low-abundance isotopomers.

We begin with the common Python practice of importing all the required packages and classes at the beginning of 
the code.

.. plot::
    :context: reset

    import matplotlib.pyplot as plt

    from mrsimulator import Simulator, Site, SpinSystem, Coupling
    from mrsimulator.method.lib import BlochDecaySpectrum
    from mrsimulator.method import SpectralDimension
    from mrsimulator import signal_processing as sp

Spin Systems
------------

The molecules in a sample of ethanol, :math:`\text{CH$_3$CH$_2$OH}`, can be formed with any of the
naturally abundant isotopes of hydrogen, carbon, and oxygen present.  Of the most abundant isotopes, 
:math:`^1\text{H}` (99.985%), :math:`^{12}\text{C}` (98.93%), and :math:`^{16}\text{O}` (99.762%), 
only :math:`^1\text{H}` is NMR active.  The most abundant NMR active isotopes of carbon and oxygen are 
:math:`^{13}\text{C}` (1.11%) and :math:`^{17}\text{O}` (0.038%).  Additionally, the 
:math:`^2\text{H}` (0.015%) isotope will be present.   For our purposes, we will ignore the effects of 
these lower abundant :math:`^{17}\text{O}` and :math:`^2\text{H}` isotopes, and focus solely on the spectra 
of the isotopomers formed from :math:`^1\text{H}`, :math:`^{12}\text{C}` , and :math:`^{13}\text{C}`.

There are three magnetically inequivalent :math:`^1\text{H}` and two magnetically inequivalent 
:math:`^{13}\text{C}` sites in ethanol.  These sites are created in the code shown below.

.. plot::
    :context: close-figs

    # All shifts in ppm
    H_CH3 = Site(isotope="1H", isotropic_chemical_shift=1.226)
    H_CH2 = Site(isotope="1H", isotropic_chemical_shift=2.61)
    H_OH = Site(isotope="1H", isotropic_chemical_shift=3.687)

    C_CH3 = Site(isotope="13C", isotropic_chemical_shift=18)
    C_CH2 = Site(isotope="13C", isotropic_chemical_shift=58)

These sites will be used, along with :ref:`coupling_documentation` objects described below, to create each of the isotopomers.

Isotopomer 1
''''''''''''

The most abundant isotopomer of ethanol consists of the :math:`^{1}\text{H}`, :math:`^{12}\text{C}`, 
and :math:`^{16}\text{O}` isotopes, as shown below.

.. figure:: ../_static/Ethanol.*
    :width: 200
    :alt: figure
    :align: center

    Most abundant isotopomer of ethanol.

Since the abundance of :math:`^{12}\text{C}` is 98.9%, the probability of this isotopomer is 
:math:`0.989 \times 0.989=0.97812`
    
Using the sites defined above, we create a list of sites present in this isotopomer.

.. plot::
    :context: close-figs
    
    iso1_sites = [H_CH3, H_CH3, H_CH3, H_CH2, H_CH2, H_OH]

Each site in the isotopomer is identified by its index in the list, which are numbered from 0 to 5.

Next we create the :ref:`coupling_documentation` objects between the sites and place the Coupling objects in
a list.

.. plot::
    :context: close-figs
    
    iso1_sites = [H_CH3, H_CH3, H_CH3, H_CH2, H_CH2, H_OH]

    # All isotropic_j shifts in ppm
    HH_coupling_1 = Coupling(site_index=[0, 3], isotropic_j=7)
    HH_coupling_2 = Coupling(site_index=[0, 4], isotropic_j=7)
    HH_coupling_3 = Coupling(site_index=[1, 3], isotropic_j=7)
    HH_coupling_4 = Coupling(site_index=[1, 4], isotropic_j=7)
    HH_coupling_5 = Coupling(site_index=[2, 3], isotropic_j=7)
    HH_coupling_6 = Coupling(site_index=[2, 4], isotropic_j=7)

    iso1_couplings = [
        HH_coupling_1,
        HH_coupling_2,
        HH_coupling_3,
        HH_coupling_4,
        HH_coupling_5,
        HH_coupling_6,
    ]

Next, we create the SpinSystem object for this isotopomer with its abundance.

.. plot::
    :context: close-figs
    
        isotopomer1 = SpinSystem(sites=iso1_sites, couplings=iso1_couplings, abundance=97.812)


Isotopomer 2
''''''''''''

Replacing the methyl carbon with a :math:`^{13}\text{C}` isotope, we get the following isotopomer pictured below (:math:`^{13}\text{C}` marked in blue)

.. figure:: ../_static/iso2.*
    :width: 200
    :alt: figure
    :align: center

    Second isotopomer of ethanol containing all :math:`^{1}\text{H}`,
    :math:`^{13}\text{C}` methyl, and :math:`^{12}\text{C}` methylene isotopes.

We now construct the spin system for this isotopomer.

.. plot::
    :context: close-figs

    iso2_sites = [H_CH3, H_CH3, H_CH3, H_CH2, H_CH2, H_OH, C_CH3]

    # Define methyl 13C - 1H couplings
    CH3_coupling_1 = Coupling(site_index=[0, 6], isotropic_j=125)
    CH3_coupling_2 = Coupling(site_index=[1, 6], isotropic_j=125)
    CH3_coupling_3 = Coupling(site_index=[2, 6], isotropic_j=125)

    # Add new couplings to existing 1H - 1H couplings
    iso2_couplings = iso1_couplings + [CH3_coupling_1, CH3_coupling_2, CH3_coupling_3]

    isotopomer2 = SpinSystem(sites=iso2_sites, couplings=iso2_couplings, abundance=1.088)

Isotopomer 3
''''''''''''

Lastly, we build the sites, couplings, and spin system for the other
isotopomer with the methylene carbon replaced with :math:`^{13}\text{C}` pictured
below (:math:`^{13}\text{C}` marked in blue)

.. figure:: ../_static/iso3.*
    :width: 200
    :alt: figure
    :align: center

    Third isotopomer of ethanol containing all :math:`^{1}\text{H}`,
    :math:`^{12}\text{C}` methyl, and :math:`^{13}\text{C}` methylene isotopes.

.. plot::
    :context: close-figs

    iso3_sites = [H_CH3, H_CH3, H_CH3, H_CH2, H_CH2, H_OH, C_CH2]

    # Define methylene 13C - 1H couplings
    CH2_coupling_1 = Coupling(site_index=[3, 6], isotropic_j=141)
    CH2_coupling_2 = Coupling(site_index=[4, 6], isotropic_j=141)

    # Add new couplings to existing 1H - 1H couplings
    iso3_couplings = iso1_couplings + [CH2_coupling_1, CH2_coupling_2]

    isotopomer3 = SpinSystem(sites=iso3_sites, couplings=iso3_couplings, abundance=1.088)




Methods
-------

Now, we define two Bloch spectrum methods for both :math:`^1\text{H}` and :math:`^{13}\text{C}`.
These methods emulate simple 1-pulse acquire experiments.

.. plot::
    :context: close-figs

    method_H = BlochDecaySpectrum(
        channels=["1H"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            SpectralDimension(
                count=16000,
                spectral_width=1.5e3,  # in Hz
                reference_offset=950,  # in Hz
                label="$^{1}$H frequency",
            )
        ],
    )

    method_C = BlochDecaySpectrum(
        channels=["13C"],
        magnetic_flux_density=9.4,  # in T
        spectral_dimensions=[
            SpectralDimension(
                count=32000,
                spectral_width=8e3,  # in Hz
                reference_offset=4e3,  # in Hz
                label="$^{13}$C frequency",
            )
        ],
    )



Simulation
----------

Now we create an instance of the simulator object, which holds a list of our three spin systems and a list of our two methods. Finally, we run the simulation.

.. plot::
    :context: close-figs

    spin_systems = [isotopomer1, isotopomer2, isotopomer3]
    methods = [method_H, method_C]
    sim = Simulator(spin_systems=spin_systems, methods=methods)
    sim.run()


Signal Processing
-----------------

Let's set up our post-simulation processing. We apply 1 Hz and 20 Hz of exponential line
broadening to the proton and carbon spectra, respectively.

.. plot::
    :context: close-figs

    # Get the simulation data
    H_data = sim.methods[0].simulation
    C_data = sim.methods[1].simulation

    # Create the signal processors
    processor_1H = sp.SignalProcessor(
        operations=[
            sp.IFFT(),
            sp.apodization.Exponential(FWHM="1 Hz"),
            sp.FFT(),
        ]
    )

    processor_13C = sp.SignalProcessor(
        operations=[
            sp.IFFT(),
            sp.apodization.Exponential(FWHM="20 Hz"),
            sp.FFT(),
        ]
    )

    # apply the signal processors
    processed_H_data = processor_1H.apply_operations(data=H_data)
    processed_C_data = processor_13C.apply_operations(data=C_data)

Plotting the Data
-----------------

Finally, we can plot the two spectra using the code below.  Additionally, we save the plot as a pdf file.

.. skip: next

.. plot::
    :context: close-figs
    :caption: :math:`^1\text{H}` and :math:`^{13}\text{C}` spectrum of ethanol. Note,
        the :math:`^{13}\text{C}` satellites seen on either side of the peaks near 1.2 ppm
        and 2.6 ppm in the :math:`^1\text{H}` spectrum.

    fig, ax = plt.subplots(
        nrows=1, ncols=2, subplot_kw={"projection": "csdm"}, figsize=[8, 3.5]
    )

    ax[0].plot(processed_H_data.real, color="black", linewidth=0.5)
    ax[0].invert_xaxis()
    ax[0].set_title("$^1$H")

    ax[1].plot(processed_C_data.real, color="black", linewidth=0.5)
    ax[1].invert_xaxis()
    ax[1].set_title("$^{13}$C")

    plt.tight_layout()
    plt.savefig("spectra.pdf")
    plt.show()


Saving your Work
----------------

If you want to save your spectrum in csdf format

.. plot::
    :context: close-figs

    processed_H_data.save("processed_H_data.csdf")
    processed_C_data.save("processed_C_data.csdf")


Saving the SpinSystems
""""""""""""""""""""""

Saving the Methods
""""""""""""""""""

Saving the full Simulation
""""""""""""""""""""""""""

.. plot::
    :include-source: False

    import os
    from os.path import isfile

    if isfile("spectra.pdf"): os.remove("spectra.pdf")
    if isfile("processed_H_data.csdf"): os.remove("processed_H_data.csdf")
    if isfile("processed_C_data.csdf"): os.remove("processed_C_data.csdf")
