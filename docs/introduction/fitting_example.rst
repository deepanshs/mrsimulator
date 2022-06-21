.. _fitting_example:

Least-Squares Fitting Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
``mrsimulator`` can interact with a variety of other Python data science 
packages.  One such package, 
`LMFIT <https://lmfit.github.io/lmfit-py/>`_, can be used to perform non-linear 
least-squares analysis of experimental NMR spectra. 

Here, we illustrate the use of the mrsimulator objects to

- create a fitting model using Simulator and SignalProcessor objects,
- use the fitting model to perform a least-squares analysis, and
- extract the fitting parameters from the model.

In this example, we apply the least-squares fitting procedure to the
:math:`^{27}\text{Al}` whole echo acquisition of :math:`\text{Al(acac)$_2$}`.

We begin by obtaining the experimental datatset measured on a 9.4 T
Bruker AVANCE IIIHD NMR spectrometer.  Bruker datasets are saved in folders
with a number as the folder name.  In this case, that folder had been 
renamed to ``Al_acac3_0``, compressed into zip archive, and uploaded to an
internet accessible server.  If you already have the folder available on your
local machine you can skip this step.

.. plot::
    :context: close-figs

    import requests
    import zipfile
    from io import BytesIO
    file_ = "https://ssnmr.org/sites/default/files/mrsimulator/Al_acac3_0.zip"
    request = requests.get(file_)
    z = zipfile.ZipFile(BytesIO(request.content))
    z.extractall("Al_acac")

.. plot::
    :context: close-figs

    import nmrglue as ng

    # initialize nmrglue converter object
    converter = ng.convert.converter()

    # read in the bruker data file
    dic, data = ng.bruker.read("Al_acac") 

    converter.from_bruker(dic,data, remove_digital_filter=False)

    # convert to CSDM format
    csdm_dataset = converter.to_csdm()

.. plot::
    :context: close-figs

    import matplotlib.pyplot as plt
    plt.figure(figsize=(5, 3))  # set the figure size
    ax = plt.subplot(projection="csdm")
    ax.plot(csdm_dataset.real)
    ax.plot(csdm_dataset.imag)
    plt.tight_layout()
    plt.grid()
    plt.show()

.. plot::
    :context: close-figs

    import numpy as np

    # set time origin to echo top
    csdm_dataset.dimensions[0].coordinates_offset = "-0.00816 s" 
    phased_dataset = csdm_dataset * np.exp(-1j* (np.pi+np.angle(csdm_dataset.max()).value))
    plt.figure(figsize=(5, 3))  # set the figure size
    ax = plt.subplot(projection="csdm")
    ax.plot(phased_dataset.real)
    ax.plot(phased_dataset.imag)
    plt.tight_layout()
    plt.grid()
    plt.show()


.. plot::
    :context: close-figs

    from mrsimulator import signal_processing as sp

    ft = sp.SignalProcessor(
        operations=[
            sp.FFT()
        ]
    )
    exp_spectrum = ft.apply_operations(data=phased_dataset)
    exp_spectrum.x[0].to("ppm", "nmr_frequency_ratio")
    plt.figure(figsize=(5, 3))  # set the figure size
    ax = plt.subplot(projection="csdm")
    ax.plot(exp_spectrum.real)
    ax.plot(exp_spectrum.imag)
    plt.tight_layout()
    ax.set_xlim(-20,20)
    plt.show()

.. plot::
    :context: close-figs

    sigma = 0.03 #guess
    exp_spectrum = exp_spectrum.real

.. plot::
    :context: close-figs

    from mrsimulator import Site, SpinSystem, Simulator
    from mrsimulator import signal_processing as sp

    site = Site(
        isotope="27Al",
        isotropic_chemical_shift=5, 
        quadrupolar = {"Cq":2.9e6, "eta":0.2},
    )
    sys = SpinSystem(sites = [site]) 

.. plot::
    :context: close-figs

    from mrsimulator.method.lib import BlochDecayCTSpectrum
    from mrsimulator.utils import get_spectral_dimensions

    spectral_dims = get_spectral_dimensions(exp_spectrum)
    MAS = BlochDecayCTSpectrum(
        channels=["27Al"],
        magnetic_flux_density=9.4,  # in T
        rotor_frequency=12500,  # in Hz
        spectral_dimensions= spectral_dims,
        experiment=exp_spectrum,  # add the measurement to the method.
    )


.. plot::
    :context: close-figs

    sim = Simulator(spin_systems=[sys], methods=[MAS])
    sim.run()

    # Post Simulation Processing
    # --------------------------
    processor = sp.SignalProcessor(
        operations=[
            sp.IFFT(),  # inverse FFT to convert frequency based spectrum to time domain.
            sp.apodization.Exponential(FWHM="50 Hz"),  # apodization of time domain signal.
            sp.FFT(),  # forward FFT to convert time domain signal to frequency spectrum.
            sp.Scale(factor=2.5e6       ),  # scale the frequency spectrum.
        ]
    )
    processed_data = processor.apply_operations(data=sim.methods[0].simulation)

    # Plot of the guess spectrum
    # --------------------------
    plt.figure(figsize=(4.25, 3.0))
    ax = plt.subplot(projection="csdm")
    ax.plot(exp_spectrum.real, "k", linewidth=1, label="Experiment")
    ax.plot(processed_data.real, "b",  linewidth=1, label="guess spectrum") #alpha=0.75,
    ax.set_xlim(-20, 20)
    plt.legend()
    plt.grid()
    plt.tight_layout()
    plt.show()


.. plot::
    :context: close-figs

    from lmfit import Minimizer
    from mrsimulator.utils import spectral_fitting as sf

    params = sf.make_LMFIT_params(sim, processor)
    params.pop("sys_0_abundance")
    print(params.pretty_print(columns=["value", "min", "max", "vary", "expr"]))

.. plot::
    :context: close-figs

    minner = Minimizer(sf.LMFIT_min_function, params, fcn_args=(sim, processor, sigma))
    result = minner.minimize()
    result

.. plot::
    :include-source: False

    import shutil

    shutil.rmtree("Al_acac")


