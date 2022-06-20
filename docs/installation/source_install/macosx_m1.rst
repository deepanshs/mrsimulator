**Arm64 Python**

The default Python installed on a M1 Mac is built both for x86_64 (Intel) and Arm64 (Apple
silicon); this will likely cause errors when running scripts. For this reason, we first need
to install Python from the ``brew`` package manager. If brew is not installed on your machine,
install brew from `homebrew <https://brew.sh>`_ before continuing.

.. note::

  To check which architecture a Python executable is built for, run the following command

  .. code-block:: bash

    file <path to python executable>

  which should output something similar to

  .. code-block:: bash

    /usr/bin/python3 (for architecture x86_64):	Mach-O 64-bit executable x86_64
    /usr/bin/python3 (for architecture arm64e):	Mach-O 64-bit executable arm64e

*Step-1* Install virtualenv

Anaconda is incompatible with Arm64 Python. We first need to install the virtualenv
package manager which is compatible with Arm64 Python.

.. code-block:: shell

  $ python3 -m pip install virtualenv

NOTE: Need to check if python or python3 by default

*Step-2* Install Arm64 Python

.. code-block:: shell

  $ brew install python@3.9

After the above command is done running, an instillation path should be printed to the console

  .. code-block:: shell

    ==> Summary
    🍺  /opt/homebrew/Cellar/python@3.9

Remember this path and replace ``<python install path>``
in the next step.

*Step-3* Create environment

  Navigate to the directory where you'd like to create a new environment folder then execute

.. code-block:: shell

  $ python3 -m virtualenv mrsimulator-env -p=<python install path>

This will create a new directory with the name ``mrsimulator-env`` which holds all of the
environments important files. To activate the environment run

  $ source ./mrsimulator-env/bin/activate

You will need to activate the environment each time the terminal is restarted.

**OpenBLAS/Accelerate and FFTW libraries**

The ``brew`` package manager is also required to install the correct the development headers
for the OpenBLAS (if applicable) and FFTW libraries.

*Step-1* Install the FFTW library using the homebrew formulae.

.. code-block:: bash

  $ brew install fftw

*Step-2* By default, the mrsimulator package links to the openblas library for BLAS
operations. Mac users may opt to choose the in-build Apple's Accelerate library. If you
opt for Apple's Accelerate library, skip to *Step-3*. If you wish to link the mrsimulator
package to the OpenBLAS library, type the following in the terminal,

.. code-block:: bash

  $ brew install openblas

*Step-3* If you choose to link the mrsimulator package to the OpenBLAS library, skip
to the next section.

*(a)* You will need to install the BLAS development header for Apple's Accelerate
library. The easiest way is to install the Xcode Command Line Tools. Note, this is a
one-time installation. If you have previously installed the Xcode Command Line Tools,
you may skip this sub-step. Type the following in the terminal,

.. code-block:: bash

  $ xcode-select --install

*(b)* The next step is to let the mrsimulator setup know your preference.
Open the ``settings.py`` file, located at the root level of the mrsimulator source
code folder, in a text editor. You should see

.. code-block:: python

  # -*- coding: utf-8 -*-
  # BLAS library
  use_openblas = True
  # mac-os only
  use_accelerate = False

To link the mrsimulator package to the Apple's Accelerate library, change the
fields to

.. code-block:: python

  # -*- coding: utf-8 -*-
  # BLAS library
  use_openblas = False
  # mac-os only
  use_accelerate = True

**Install a C/C++ compiler**

The C-compiler installs with the Xcode Command Line Tools. No further action is
required.
