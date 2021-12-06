.. _method_documentation:

=================
The Method Object
=================

Here is where the documentation and example for how to use the :ref:`method_api` object will
eventually go. For now this page is just a placeholder.

Need to document Named Methods, simulation & experiment, SpectralDimensions, SpectralEvent, MixingEvent, TransitionQuery, MixingQuery,

A :ref:`method_api` object is a collection of attributes that describe an NMR method.
In ``mrsimulator``, all methods are described through five keywords -

In ``mrsimulator`` we have two classes of methods. The first are **named methods** and the second
are **custom methods**

Named Methods
'''''''''''''

Named methods are methods with some pre-defined attributes

Custom Methods
''''''''''''''

Custom methods can contain whatever you like

Reference Tables
----------------

.. cssclass:: table-bordered table-striped centered
.. _table_method:
.. list-table:: The attributes of a Method, Method1D, and Method2D object
  :widths: 15 20 65
  :header-rows: 1

  * - Attribute Name
    - Type
    - Description

  * - channels
    - ``List``
    - A *required* list of isotopes given as strings over which the given method applies.
      For example, ``["1H"]``.

  * - magnetic_flux_density
    - ``float``
    - An *optional* float describing the macroscopic magnetic flux density of the applied
      external magnetic filed in Tesla. For example, ``18.8`` Tesla. The default is ``9.4`` tesla.

  * - spectral_dimensions
    - ``List`` of :ref:`spectral_dim_api` objects
    - A list of spectral dimensions. The coordinates along each spectral dimension is
      described with the keywords, *count* (:math:`N`), *spectral_width*
      (:math:`\nu_\text{sw}`), and *reference_offset* (:math:`\nu_0`). The
      coordinates are evaluated as,

      .. math::
        \left([0, 1, 2, ... N-1] - \frac{T}{2}\right) \frac{\nu_\text{sw}}{N} + \nu_0

      where :math:`T=N` when :math:`N` is even else :math:`T=N-1`.



.. cssclass:: table-bordered

.. list-table::
  :widths: 25 75
  :header-rows: 1

  * - Keywords
    - Description
  * - channels
    - A list of isotope symbols over which the given method applies.
  * - magnetic_flux_density
    - The macroscopic magnetic flux density of the applied external magnetic field.
  * - rotor_angle
    - The angle between the sample rotation axis and the applied external magnetic field.
  * - rotor_frequency
    - The sample rotation frequency.
  * - spectral_dimensions
    - A list of spectral dimensions. The coordinates along each spectral dimension is
      described with the keywords, *count* (:math:`N`), *spectral_width*
      (:math:`\nu_\text{sw}`), and *reference_offset* (:math:`\nu_0`). The
      coordinates are evaluated as,

      .. math::
        \left([0, 1, 2, ... N-1] - \frac{T}{2}\right) \frac{\nu_\text{sw}}{N} + \nu_0

      where :math:`T=N` when :math:`N` is even else :math:`T=N-1`.
