
**OpenBLAS and FFTW libraries**

On Linux, the package manager for your distribution is usually the easiest route to
ensure you have the prerequisites to building the **mrsimulator** library. To build from
source, you will need the OpenBLAS and FFTW development headers for your Linux
distribution. Type the following command in the terminal, based on your Linux
distribution.

*For (Debian/Ubuntu):*

.. code-block:: bash

  $ sudo apt-get install libopenblas-dev libfftw3-dev

*For (Fedora/RHEL):*

.. code-block:: bash

  $ sudo yum install openblas-devel fftw-devel

**Install a C/C++ compiler**

The C-compiler comes with your Linux distribution. No further action is
required.
