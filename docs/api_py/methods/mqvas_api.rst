.. currentmodule:: mrsimulator.methods

.. _mqvas_ref:

Multi-quantum variable-angle spinning
-------------------------------------

The following classes are used in simulating multi-quantum variable-angle spinning
spectrum correlating the frequencies from the symmetric multiple-quantum transition to
the central transition frequencies. The :math:`p` and :math:`d` pathways for the MQVAS
methods are

.. math::
    p: &0 \rightarrow M \rightarrow -1 \\
    d: &0 \rightarrow 0 \rightarrow 0

where :math:`M` is the multiple-quantum number. The value of :math:`M` depends on the
spin quantum number, :math:`I`, and is listed in :numref:`tb_mqvas`.

**Affine mapping**

The resulting spectrum is sheared and scaled, such that the frequencies along the indirect
dimension are given as

.. math::
    \langle \Omega\rangle_\text{MQ-VAS} =  \frac{1}{1+\kappa}\Omega_{m, -m} +
            \frac{\kappa}{1+\kappa}\Omega_{1/2, -1/2}.

Here, :math:`\langle \Omega\rangle_\text{MQ-VAS}` is the average frequency along the indirect
dimension, :math:`\Omega_{m, -m}` and :math:`\Omega_{1/2, -1/2}` are the frequency
contributions from the :math:`|m\rangle \rightarrow |-m\rangle` symmetric multiple-quantum
transition and the central transition, respectively, and :math:`\kappa` is the shear factor.
The values of the shear factor for various transitions are listed in :numref:`tb_mqvas`.

.. cssclass:: table-bordered

.. _tb_mqvas:
.. list-table:: The table lists the multi-quantum transition associated with the spin :math:`I`,
        and the corresponding shear factor, :math:`\kappa`, used in affine mapping of the MQ-VAS
        methods.
    :widths: 20 40 20 20
    :header-rows: 1

    * - Spin
      - Symmetric multi-quantum transition
      - :math:`M`
      - :math:`\kappa`

    * - 3/2
      - :math:`\left(\frac{3}{2} \rightarrow -\frac{3}{2}\right)`
      - :math:`-3`
      - 21/27

    * - 5/2
      - :math:`\left(-\frac{3}{2} \rightarrow \frac{3}{2}\right)`
      - :math:`3`
      - 114/72

    * - 5/2
      - :math:`\left(\frac{5}{2} \rightarrow -\frac{5}{2}\right)`
      - :math:`-5`
      - 150/72

    * - 7/2
      - :math:`\left(-\frac{3}{2} \rightarrow \frac{3}{2}\right)`
      - :math:`3`
      - 303/135

    * - 7/2
      - :math:`\left(-\frac{5}{2} \rightarrow \frac{5}{2}\right)`
      - :math:`5`
      - 165/135

    * - 7/2
      - :math:`\left(\frac{7}{2} \rightarrow -\frac{7}{2}\right)`
      - :math:`-7`
      - 483/135

    * - 9/2
      - :math:`\left(-\frac{3}{2} \rightarrow \frac{3}{2}\right)`
      - :math:`3`
      - 546/216

    * - 9/2
      - :math:`\left(-\frac{5}{2} \rightarrow \frac{5}{2}\right)`
      - :math:`5`
      - 570/216

    * - 9/2
      - :math:`\left(-\frac{7}{2} \rightarrow \frac{7}{2}\right)`
      - :math:`5`
      - 84/216

.. _threeQ_vas_ref:

Triple-quantum variable-angle spinning method
'''''''''''''''''''''''''''''''''''''''''''''

.. autoclass:: ThreeQ_VAS

.. _fiveQ_vas_ref:

Five-quantum variable-angle spinning method
'''''''''''''''''''''''''''''''''''''''''''

.. autoclass:: FiveQ_VAS

.. _sevenQ_vas_ref:

Seven-quantum variable-angle spinning method
''''''''''''''''''''''''''''''''''''''''''''

.. autoclass:: SevenQ_VAS
