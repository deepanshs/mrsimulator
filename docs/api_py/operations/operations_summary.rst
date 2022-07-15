.. _operations_api:

Operations
==========

Generic operations
------------------

.. currentmodule:: mrsimulator.signal_processor

Import the module as

.. code-block:: python

    from mrsimulator import signal_processor as sp

.. rubric:: Operation Summary

The following list of operations applies to **all dependent variables** within the
CSDM object.

.. autosummary::
    :nosignatures:

      ~Scale
      ~Linear
      ~IFFT
      ~FFT

Baseline
--------

.. currentmodule:: mrsimulator.signal_processor.baseline

Access the sub-module as ``sp.baseline``

.. rubric:: Operation Summary

The following list of operations applies to **selected dependent variables** within
the CSDM object.

.. autosummary::
    :nosignatures:

      ~Polynomial
      ~ConstantOffset

.. seealso::

    :ref:`signal_processor_documentation` for more details.

    :ref:`signal_processor_examples` for notebooks using these operations.


Apodization
-----------

.. currentmodule:: mrsimulator.signal_processor.apodization

Access the sub-module as ``sp.apodization``

.. rubric:: Operation Summary

The following list of operations applies to **selected dependent variables** within
the CSDM object.

.. autosummary::
    :nosignatures:

      ~Gaussian
      ~Exponential
      ~SkewedGaussian
      ~TopHat
      ~Mask

.. seealso::

    :ref:`signal_processor_documentation` for more details.

    :ref:`signal_processor_examples` for notebooks using these operations.


Affine Transformation
---------------------

.. currentmodule:: mrsimulator.signal_processor.affine

Access the sub-module as ``sp.affine``

.. rubric:: Operation Summary

The following list of operations applies to **selected dependent variables** within
the CSDM object.

.. autosummary::
    :nosignatures:

      ~Shear
      ~Scale

.. seealso::

    :ref:`signal_processor_documentation` for more details.
