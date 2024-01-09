
**Install conda**

Skip this step if you already have miniconda or anaconda for python>=|PY_VERSION| installed on
your system.
Download the latest version of conda on your operating system from either
`miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ or
`Anaconda <https://www.anaconda.com/products/individual/>`_ websites. Make sure you
download conda for Python 3. Double click the downloaded .exe file and follow the
installation steps.

**OpenBLAS and FFTW libraries**

Launch the ``Anaconda prompt`` (it should be located under the start menu). Within the
anaconda prompt, type the following to install the package dependencies.

.. code-block:: bash

    $ conda install -c conda-forge openblas fftw

**Install a C/C++ compiler**

Because the core of the **mrsimulator** package is written in C, you will require a
C-compiler to build and install the package. Download and install the Microsoft
Visual C++ compiler from
`Build Tools for Visual Studio <https://visualstudio.microsoft.com/visual-cpp-build-tools/>`_.
