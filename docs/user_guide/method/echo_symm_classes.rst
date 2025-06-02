.. note::

    The well-known Hahn-echo can occur whenever the :math:`p_I` values of
    transitions in a transition pathway change sign.  This is because the
    changing sign of :math:`p_I` leads to a sign change for every
    :math:`p_I`-dependent transition frequency contribution. Thus, a Hahn
    echo forms whenever

    .. math::
        \overline{\text{p}_I} = \frac{1}{t} \int_0^t \text{p}_I(t') \, dt' = 0,

    assuming a frequency contribution's spatial symmetry function, :math:`{\Xi}`,
    remains constant during this period.  As seen in the table in the
    :ref:`frequency_contribution_table` table, sign changes in other symmetry
    functions can also lead to corresponding sign changes for dependent
    frequency contributions.  Thus, a problem with showing only the :math:`p_I`
    symmetry pathway for an NMR method is that it does not explain the formation
    of other classes of echoes that result when other symmetry functions change
    sign in a transition pathway.  To fully understand when and which frequency
    contributions refocus into echoes, we must follow *all* relevant spatial,
    transition, or spatial-transition product symmetries through an NMR
    experiment.   Thus, we generally classify echoes that refocus during a time
    interval as a *transition symmetry echo* (at constant :math:`{\Xi}_k`) when

    .. math::
        \overline{{\xi}_k} = \frac{1}{t} \int_0^t {\xi}_k(t') \, dt' = 0,

    and as a *spatial symmetry echo* (at constant :math:`{\xi}_k`) when

    .. math::
        \overline{{\Xi}_k} = \frac{1}{t} \int_0^t {\Xi}_k(t') \, dt' = 0,

    and as a *spatial-transition symmetry product* echo when

    .. math::

        \overline{{\Xi}_k {\xi}_k} = \frac{1}{t} \int_0^t {\Xi}_k(t') \, {\xi}_k(t')  \, dt' = 0.

    Within the class of transition echoes, we find subclasses such as
    :math:`\text{p}` echoes, which include the Hahn echo and the stimulated
    echo; :math:`\text{d}` echoes, which include the solid echo and Solomon
    echoes,  :math:`\text{c}_4` echoes, used in MQ-MAS and Satellite-Transition
    Magic-Angle Spinning (ST-MAS); :math:`\text{c}_2` echoes, used in
    Correlation Of Anisotropies Separated Through Echo Refocusing (COASTER); and
    :math:`\text{c}_0` echoes, used in Multiple-Quantum DOuble Rotation
    (MQ-DOR).

    Within the class of spatial echoes, we find subclasses such as :math:`\mathbb{D}`
    rotary echoes, which occur during sample rotation, and :math:`\mathbb{D}_0` and
    :math:`\mathbb{G}_0` echoes, which are designed to occur simultaneously during the
    Dynamic-Angle Spinning (DAS) experiment.
