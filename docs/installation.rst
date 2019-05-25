

.. _shielding_tensor_api:

============
Installation
============

The ``mrsimulator`` package requires the 'fftw3' package to function properly.
Install the fftw3 package using

.. doctest::
    :skipif: None is None

    >>> conda install -c eumetsat fftw3

Read more about
`fftw3 installation <https://anaconda.org/eumetsat/fftw3>`_.


Before installing the package, install the required dependency packages using

.. doctest::
    :skipif: None is None

    >>> pip install scipy numpy astropy mkl mkl-include

To install the ``mrsimulator`` package, first download package and run the
following in the terminal

.. doctest::
    :skipif: None is None

    >>> python setup.py install
