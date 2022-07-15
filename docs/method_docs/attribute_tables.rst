.. _method_attribute_tables:

===================
Attribute Summaries
===================

.. cssclass:: table-bordered table-striped centered
.. _table_method:
.. list-table:: The attributes of a Method object
  :widths: 20 15 65
  :header-rows: 1

  * - Attribute Name
    - Type
    - Description

  * - channels
    - ``list``
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
    - ``list``
    - A list of :ref:`spectral_dim_api` objects describing the spectral dimensions for the method.

  * - affine_matrix
    - ``list`` or ``np.ndarray``
    - An *optional* (``n`` x ``n``) affine transformation matrix represented by a numpy array where ``n`` is
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
    - An *optional* float representing the origin offset, or Larmor frequency, along the
      spectroscopic dimension in units of Hz. The default value is ``None`` and the origin offset
      is set to the Larmor frequency of isotope from the :attr:`~mrsimulator.Method.channels`
      attribute of the method containing the spectral dimension.

  * - events
    - ``List``
    - An *optional* list of :ref:`event_api` used to emulate an experiment.
      The default value is an empty list.


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
      ``None`` and takes the global magnetic flux density defined by the method's
      :attr:`~mrsimulator.Method.magnetic_flux_density` attribute.

  * - rotor_angle
    - ``float``
    - An *optional* float describing the angle between the sample rotation axis and the external
      magnetic field in radians. The default is ``None`` and takes the global rotor angle defined
      by the method's :attr:`~mrsimulator.Method.rotor_angle` attribue.

  * - rotor_frequency
    - ``float``
    - An *optional* float describing the sample rotation frequency in Hz. For example, ``2000`` Hz.
      The default value is ``None`` and takes the global rotor frequency defined by the method's
      :attr:`~mrsimulator.Method.rotor_frequency` attribute.

  * - freq_contrib
    - ``List``
    - An *optional* list of :ref:`freq_contrib_api` (list of allowed strings) selecting which
      contributions to include when calculating a transition frequency. For example,
      ``["Shielding1_0", "Shielding1_2"]``.  When unspecified, it defaults to a list holding 
      all frequency enumerations.

  * - transition_queries
    - ``list``
    - An *optional* ``list`` of :ref:`query_api` objects, or their ``dict`` representations
      selecting transitions active during the event. The default is one **TransitionQuery** 
      with *P=[0]* on ``ch1`` and ``None`` on all other channels


.. cssclass:: table-bordered table-striped centered
.. _table_mixing_event:
.. list-table:: The attributes of a MixingEvent object
  :widths: 20 15 65
  :header-rows: 1

  * - Attribute Name
    - Type
    - Description

  * - query
    - ``dict`` or :py:class:`~mrsimulator.method.MixingQuery`
    - A :py:class:`~mrsimulator.method.MixingQuery` object, or its ``dict`` representation,
      determines the complex amplitude of mixing between transitions in adjacent spectral or delay
      events.

..   - The coordinates along each spectral dimension are
..       described with the keywords,``count``(:math:`N`), ``spectral_width``
..       (:math:`\nu_\text{sw}`), and``reference_offset``(:math:`\nu_0`). The
..       coordinates are evaluated as,
..
..       .. math
..         \left([0, 1, 2, ... N-1] - \frac{T}{2}\right) \frac{\nu_\text{sw}}{N} + \nu_0
..
..       where :math:`T=N` when :math:`N` is even else :math:`T=N-1`.
