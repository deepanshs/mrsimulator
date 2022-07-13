.. _origin_and_reference_offset_documentation:

===========================
Origin and Reference Offset
===========================

:py:meth:`~mrsimulator.method.spectral_dimension.SpectralDimension` has additional
attributes that have already been discussed in earlier sections of the documentation.
Notably, ``origin_offset`` and ``reference_offset`` are important for converting
the frequency coordinate into a dimensionless frequency ratio coordinate. For
spectra where all the spectral dimensions are associated with single-quantum
transitions on a single isotope, the convention for defining ``origin_offset``
and ``reference_offset`` is well established;
the ``origin_offset``, :math:`o_k`, is interpreted as the NMR spectrometer
frequency and  the ``reference_offset``, :math:`b_k`, as the reference
frequency. Given the frequency coordinate, :math:`{X}`, the corresponding
dimensionless-frequency ratio follows,

.. math::
    :label: chemicalShiftDef

    {X}^\text{ratio} = \displaystyle \frac{{X}}{o_k - b_k}.

In the case of multiple quantum dimensions, however, there appear
to be no formal conventions for defining ``origin_offset`` and ``reference_offset``.
