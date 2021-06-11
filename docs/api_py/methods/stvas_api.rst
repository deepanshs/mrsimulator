.. currentmodule:: mrsimulator.methods

.. _stvas_ref:

Satellite-transition variable-angle spinning (ST-VAS)
"""""""""""""""""""""""""""""""""""""""""""""""""""""

The following classes are used in simulating satellite-transition variable-angle spinning
spectrum correlating the frequencies from the satellite transitions to the central transition
frequencies. The :math:`p` and :math:`d` pathways for the ST-VAS methods are

.. math::
    p: &0 \rightarrow -1 \rightarrow -1 \\
    d: &0 \rightarrow \pm d_0 \rightarrow 0

where :math:`d_0 = m_f^2 -  m_i^2` for transition :math:`|m_i\rangle \rightarrow |m_f\rangle`.
The value of :math:`n` depends on the spin quantum number, :math:`I`, and is listed in
:numref:`tb_stvas`.

**Affine mapping**

The resulting spectrum is sheared and scaled, such that the frequencies along indirect
dimension are given as

.. math::
    \langle \Omega\rangle_\text{ST-VAS} = \frac{1}{1+\kappa}\Omega_{m, m-1} +
            \frac{\kappa}{1+\kappa}\Omega_{1/2, -1/2}.

Here, :math:`\langle \Omega\rangle_\text{ST-VAS}` is the average frequency along the indirect
dimension, :math:`\Omega_{m, m-1}` and :math:`\Omega_{1/2, -1/2}` are the frequency
contributions from the :math:`|m\rangle \rightarrow |m-1\rangle` satellite transition and the
central transition, respectively, and :math:`\kappa` is the shear factor. The values of the
shear factor for various satellite transitions are listed in :numref:`tb_stvas`.


.. cssclass:: table-bordered
.. _tb_stvas:
.. list-table:: The table lists the satellite transitions associated with the spin :math:`I`,
        and the corresponding shear factor, :math:`\kappa`, used in affine mapping of the ST-VAS
        methods.
    :widths: 20 40 20 20
    :header-rows: 1

    * - Spin
      - Satellite transitions
      - :math:`d_0`
      - :math:`\kappa`

    * - 3/2
      - :math:`\left(\frac{3}{2} \rightarrow \frac{1}{2}\right)`, :math:`\left(-\frac{1}{2} \rightarrow -\frac{3}{2}\right)`
      - :math:`2`
      - 24/27

    * - 5/2
      - :math:`\left(-\frac{3}{2} \rightarrow -\frac{1}{2}\right)`, :math:`\left(\frac{1}{2} \rightarrow \frac{3}{2}\right)`
      - :math:`2`
      - 21/72

    * - 5/2
      - :math:`\left(\frac{5}{2} \rightarrow \frac{3}{2}\right)`, :math:`\left(-\frac{3}{2} \rightarrow -\frac{5}{2}\right)`
      - :math:`4`
      - 132/72

    * - 7/2
      - :math:`\left(-\frac{3}{2} \rightarrow -\frac{1}{2}\right)`, :math:`\left(\frac{1}{2} \rightarrow \frac{3}{2}\right)`
      - :math:`2`
      - 84/135

    * - 7/2
      - :math:`\left(-\frac{5}{2} \rightarrow -\frac{3}{2}\right)`, :math:`\left(\frac{3}{2} \rightarrow \frac{5}{2}\right)`
      - :math:`4`
      - 69/135

    * - 9/2
      - :math:`\left(-\frac{3}{2} \rightarrow -\frac{1}{2}\right)`, :math:`\left(\frac{1}{2} \rightarrow \frac{3}{2}\right)`
      - :math:`2`
      - 165/216

    * - 9/2
      - :math:`\left(-\frac{5}{2} \rightarrow -\frac{3}{2}\right)`, :math:`\left(\frac{3}{2} \rightarrow \frac{5}{2}\right)`
      - :math:`4`
      - 12/216

.. _st1_vas_ref:

Inner satellite variable-angle spinning method
''''''''''''''''''''''''''''''''''''''''''''''

.. autoclass:: ST1_VAS
    :show-inheritance:
    :members:
    :inherited-members: BaseModel

.. _st2_vas_ref:

Second to inner satellite variable-angle spinning method
''''''''''''''''''''''''''''''''''''''''''''''''''''''''

.. autoclass:: ST2_VAS
    :show-inheritance:
    :members:
    :inherited-members: BaseModel
