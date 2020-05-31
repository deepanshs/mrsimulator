
**Step-1** Install conda. (Skip this step if you already have miniconda or anaconda
for python>=3.6 installed on your system.)

Download the latest version of conda on your operating system from
`miniconda <https://docs.conda.io/en/latest/miniconda.html>`_. Make sure you download
conda for python 3. Double click the downloaded .exe file and follow the installation
steps.

**Step-2** Install dependencies using conda.

Launch the ``Anaconda prompt`` (it should be located under the start menu). Within the
anaconda prompt, type the following to install the package dependencies.

.. code-block:: bash

    $ conda install -c conda-forge openblas fftw>=3.3.0

**Step-3** Install a C/C++ compiler.

Because the core of the mrsimulator package is written in C, you will require a
C-compiler to build and install the package. Download and install the Microsoft
Visual C++ compiler from
`Build Tools for Visual Studio 2019 <https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2019>`_.
