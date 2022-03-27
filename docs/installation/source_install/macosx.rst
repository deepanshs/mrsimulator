
**OpenBLAS/Accelerate and FFTW libraries**

You will require the ``brew`` package manager to install the development headers for the
OpenBLAS (if applicable) and FFTW libraries. Read more on installing brew from
`homebrew <https://brew.sh>`_.

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
