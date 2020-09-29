.. _methods_summary_api:

Methods
=======

.. currentmodule:: mrsimulator.methods

The following are the list of methods currently supported by ``mrsimulator`` as a part
of the ``mrsimulator.methods`` module. To import a method, for example the
`BlochDecaySpectrum`, used

.. doctest::

    >>> from mrsimulator.methods import BlochDecaySpectrum

All methods categorize into two groups, generic and specialized methods. A generic
method is general and bases itself on the number of spectral dimensions. At present,
there are two generic methods. All specialized methods are derived from their respective
generic method objects. The purpose of the specialized methods is to facilitate user ease
when setting up some commonly used methods, such as MQMAS, STMAS, PASS, MAT, etc.


Generic methods
---------------

.. autosummary::
    ~Method1D
    ~Method2D

Specialized methods
-------------------

Specialized one-dimensional methods
'''''''''''''''''''''''''''''''''''

.. autosummary::
    ~BlochDecaySpectrum
    ~BlochDecayCentralTransitionSpectrum


Specialized two-dimensional methods
'''''''''''''''''''''''''''''''''''

.. _mqvas_ref:

Multi-quantum variable-angle spinning (MQVAS)
"""""""""""""""""""""""""""""""""""""""""""""

The following classes are used in simulating multi-quantum variable-angle spinning
spectrum correlating the frequencies from the symmetric multiple-quantum transition to
the central transition frequencies. The :math:`p` and :math:`d` pathways for the MQVAS
methods are

.. math::
      \begin{align}
          p: &0 \rightarrow M \rightarrow -1 \\
          d: &0 \rightarrow 0 \rightarrow 0
      \end{align},

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


.. autosummary::
    ~ThreeQ_VAS
    ~FiveQ_VAS
    ~SevenQ_VAS


.. _stvas_ref:

Satellite-transition variable-angle spinning (ST-VAS)
"""""""""""""""""""""""""""""""""""""""""""""""""""""

The following classes are used in simulating satellite-transition variable-angle spinning
spectrum correlating the frequencies from the satellite transitions to the central transition
frequencies. The :math:`p` and :math:`d` pathways for the ST-VAS methods are

.. math::
      \begin{align}
          p: &0 \rightarrow -1 \rightarrow -1 \\
          d: &0 \rightarrow \pm d_0 \rightarrow 0
      \end{align},

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
shear factor for various satellite transitions are listed in :numref:`tb_stvas`..


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


.. autosummary::
    ~ST1_VAS
    ~ST2_VAS
