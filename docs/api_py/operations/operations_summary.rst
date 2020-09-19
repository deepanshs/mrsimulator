.. _operations_api:

Operations
==========

Generic operations
------------------

.. currentmodule:: mrsimulator.signal_processing

Import the module as

.. doctest::

    >>> import mrsimulator.signal_processing as sp

.. rubric:: Operation Summary

The following list of operations applies to **all dependent variables** within the
CSDM object.

.. autosummary::
    :nosignatures:

      ~Scale
      ~IFFT
      ~FFT

Apodization
-----------

.. currentmodule:: mrsimulator.signal_processing.apodization

Import the module as

.. doctest::

    >>> import mrsimulator.signal_processing.apodization as apo

.. rubric:: Operation Summary

The following list of operations applies to **selected dependent variables** within
the CSDM object.

.. autosummary::
    :nosignatures:

      ~Gaussian
      ~Exponential

.. seealso::

    :ref:`signal_processing` for a details.


Affine Transformation
---------------------

.. currentmodule:: mrsimulator.signal_processing.affine

Import the module as

.. doctest::

    >>> import mrsimulator.signal_processing.affine as af

.. rubric:: Operation Summary

The following list of operations applies to **selected dependent variables** within
the CSDM object.

.. autosummary::
    :nosignatures:

      ~Shear
      ~Scale

.. seealso::

    :ref:`signal_processing` for a details.
