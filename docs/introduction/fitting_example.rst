.. _fitting_example:

Least-Squares Fitting Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
``mrsimulator`` can also operate with other Python packages, such as 
`LMFIT <https://lmfit.github.io/lmfit-py/>`_,  to perform non-linear 
least-squares analysis of experimental NMR spectra. 

Here, we illustrate the use of the mrsimulator objects to

- create a fitting model using Simulator and SignalProcessor objects,
- use the fitting model to perform a least-squares analysis, and
- extract the fitting parameters from the model.

In this example, we apply the least-squares fitting procedure to the
:math:`^{27}\text{Al}` whole echo acquisition of :math:`\text{Al(acac)$_2$}`.

We begin by obtaining the experimental datatset measured on a 9.4 T
Bruker AVANCE IIIHD NMR spectrometer.  Bruker datasets are saved in folders
with some number as the folder name.  In this case, that folder had been 
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
    import shutil

    # initialize nmrglue converter object
    converter = ng.convert.converter()

    # load data into the converter
    dic, data = ng.bruker.read("Al_acac") # read in the bruker data file

    shutil.rmtree("Al_acac")
    converter.from_bruker(dic,data, remove_digital_filter=True)

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
    plt.show()

.. plot::
    :context: close-figs

    import numpy as np

    # set time origin to echo top
    csdm_dataset.dimensions[0].coordinates_offset = "-0.008 s" 
    phased_dataset = csdm_dataset * np.exp(-1j* np.angle(csdm_dataset.max()).value)*np.exp(-1j*np.pi)
    plt.figure(figsize=(5, 3))  # set the figure size
    ax = plt.subplot(projection="csdm")
    ax.plot(phased_dataset.real)
    ax.plot(phased_dataset.imag)
    plt.tight_layout()
    plt.show()


.. plot::
    :context: close-figs

    from mrsimulator import signal_processing as sp

    ft = sp.SignalProcessor(
        operations=[
            sp.FFT()
        ]
    )
    freq_data = ft.apply_operations(data=phased_dataset)
    freq_data.x[0].to("ppm", "nmr_frequency_ratio")
    plt.figure(figsize=(5, 3))  # set the figure size
    ax = plt.subplot(projection="csdm")
    ax.set_xlim(-20,20)
    ax.plot(freq_data.real)
    ax.plot(freq_data.imag)
    plt.tight_layout()
    plt.show()


