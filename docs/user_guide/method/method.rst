.. _method_documentation:

======
Method
======

Here is where the documentation and example for how to use the :ref:`method_api` object will
eventually go. For now this page is just a placeholder.

Need to document Named Methods, simulation & experiment, SpectralDimensions, SpectralEvent, MixingEvent, TransitionQuery, MixingQuery,

Include note about global attributes and broadcasting to events

A :ref:`method_api` object is a collection of attributes that describe an NMR method.
In ``mrsimulator``, all methods are described through five keywords -

Reference Tables
----------------

.. cssclass:: table-bordered table-striped centered
.. _table_method:
.. list-table:: The attributes of a Method, Method1D, and Method2D object
  :widths: 20 15 65
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
      external magnetic field in tesla. For example, ``18.8`` tesla. The default value is
      ``9.4`` tesla.

  * - rotor_frequency
    - ``float``
    - An *optional* float describing the sample rotation frequency in Hz. For example, ``2000`` Hz.
      The default value is ``0`` Hz.

  * - rotor_angle
    - ``float``
    - An *optional* float describing the angle between the sample rotation axis and the external
      magnetic field in radians. The default value is the magic angle,
      ``54.735 * 3.14159 / 180 = 0.955305`` radians.

  * - spectral_dimensions
    - ``List``
    - A list of :ref:`spectral_dim_api` objects describing the spectral dimensions for the method.
      For more discussion on spectral dimension objects, see the :ref:`spectral_dim_documentation`
      documentation ((NOT WRITTEN YET)).

  * - affine_matrix
    - ``np.ndarray``
    - A (``n`` x ``n``) affine transformation matrix represented by a numpy array where ``n`` is
      the number of spectral dimensions. If provided, the transformation is applied after running
      a simulation. The default value is ``None`` and no transformation is applied.

  * - simulation
    - CSDM object
    - A CSDM object representing the spectrum simulated by the method. By default, the value is
      ``None``. A value is assigned to this attribute when you run the
      simulation using the :py:meth:`~mrsimulator.Simulator.run` method.

  * - experiment
    - CSDM object
    - An *optional* CSDM object holding an experimental measurement of the method. The default
      value is ``None``


.. cssclass:: table-bordered table-striped centered
.. _table_spectral_dim:
.. list-table:: The attributes of a SpectralDimension object
  :widths: 20 15 65
  :header-rows: 1

  * - Attribute Name
    - Type
    - Description

  * - count
    - ``int``
    - An *optional* integer representing the number of points, :math:`N`, along the spectroscopic
      dimension. For example, ``4096``. The default value is ``1024``.

  * - spectral_width
    - ``float``
    - An *optional* float representing the width, :math:`\Delta x`, of the spectroscopic dimension
      in Hz. For example, ``10e3`` for 10 kHz. The default value is ``25000`` Hz.

  * - reference_offset
    - ``float``
    - An *optional* float representing the reference offset, :math:`x_0`, of the spectroscopic
      dimension in Hz. For example, ``-8000`` Hz. The default value is ``0``.

  * - origin_offset
    - ``float``
    - An optional float representing the origin offset, or Larmor frequency, along the
      spectroscopic dimension in units of Hz. The default value is ``None`` and the origin offset
      is set to the Larmor frequency of isotope from the :attr:`~mrsimulator.Method.channels`
      attribute of the method containing the spectral dimension.

  * - events
    - ``List``
    - An *optional* list of :ref:`event_api` objects used to emulate an experiment. For more
      discussion on event objects, see the :ref:`event_documentation` ((NOT WRITTEN YET)).
      The default value is a list with a single **SpectralEvent** with a symmetry_query of
      P=[-1]


.. cssclass:: table-bordered table-striped centered
.. _table_spectral_event:
.. list-table:: The attributes of a SpectralEvent object
  :widths: 20 15 65
  :header-rows: 1

  * - Attribute Name
    - Type
    - Description

  * - magnetic_flux_density
    - ``float``
    - An *optional* float describing the macroscopic magnetic flux density of the applied
      external magnetic field in tesla. For example, ``18.8`` tesla. The default value is
      ``None`` and takes the global magnetic flux density defined by
      :attr:`~mrsimulator.Method.magnetic_flux_density`.

  * - rotor_angle
    - ``float``
    - An *optional* float describing the angle between the sample rotation axis and the external
      magnetic field in radians. The default is ``None`` and takes the global rotor angle defined
      by :attr:`~mrsimulator.Method.rotor_angle`.

  * - rotor_frequency
    - ``float``
    - An *optional* float describing the sample rotation frequency in Hz. For example, ``2000`` Hz.
      The default value is ``None`` and takes the global rotor frequency defined by
      :attr:`~mrsimulator.Method.rotor_frequency`.

  * - freq_contrib
    - ``List``
    - An *optional* list of :ref:`freq_contrib_api` ((object?)) selecting which frequency
      contributions to include when calculating the spectrum. For example,
      ``["Shielding1_0", "Shielding1_2"]``. By default, the list is all frequency enumerations and
      all frequency contributions are calculated.

  * - transition_query
    - ``dict`` or :ref:`transition_api`
    - An *optional* ``dict`` or :ref:`transition_api` selecting transitions active
      during the event. Only these selected transitions will contribute to the net frequency.


.. cssclass:: table-bordered table-striped centered
.. _table_mixing_event:
.. list-table:: The attributes of a MixingEvent object
  :widths: 20 15 65
  :header-rows: 1

  * - Attribute Name
    - Type
    - Description

  * - mixing_query
    - ``dict``
    - A mixing_query object ((need to decide description of object))

..   - The coordinates along each spectral dimension is
..       described with the keywords, *count* (:math:`N`), *spectral_width*
..       (:math:`\nu_\text{sw}`), and *reference_offset* (:math:`\nu_0`). The
..       coordinates are evaluated as,
..
..       .. math
..         \left([0, 1, 2, ... N-1] - \frac{T}{2}\right) \frac{\nu_\text{sw}}{N} + \nu_0
..
..       where :math:`T=N` when :math:`N` is even else :math:`T=N-1`.
