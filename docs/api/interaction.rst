

====
Core
====

.. _haeberlen_api:

--------------------
Haeberlen Convension
--------------------

.. currentmodule:: mrsimulator.core

.. autoclass:: HaeberlenConvension
    :show-inheritance:

    .. rubric:: Attributes Summary

    .. autosummary::
        ~HaeberlenConvension.iso
        ~HaeberlenConvension.zeta
        ~HaeberlenConvension.eta

    .. rubric:: Methods Summary

    .. autosummary::
        ~HaeberlenConvension.principal_values
        ~HaeberlenConvension.dumps

    .. rubric:: Attributes Documentation
    .. autoattribute:: iso
    .. autoattribute:: zeta
    .. autoattribute:: eta

    .. rubric:: Method Documentation
    .. automethod:: principal_values
    .. automethod:: dumps

.. _principal_api:

--------------------
Principal Components
--------------------

.. currentmodule:: mrsimulator.core

.. autoclass:: PrincipalComponents
    :show-inheritance:

    .. rubric:: Attributes Summary

    .. autosummary::
        ~PrincipalComponents.xx
        ~PrincipalComponents.yy
        ~PrincipalComponents.zz

    .. rubric:: Methods Summary

    .. autosummary::
        ~PrincipalComponents.haeberlen_values
        ~PrincipalComponents.dumps

    .. rubric:: Attributes Documentation
    .. autoattribute:: xx
    .. autoattribute:: yy
    .. autoattribute:: zz

    .. rubric:: Method Documentation
    .. automethod:: haeberlen_values
    .. automethod:: dumps


.. _euler_angle_api:

------------
Euler Angles
------------

.. currentmodule:: mrsimulator.core

.. autoclass:: EulerAngles
    :show-inheritance:

    .. rubric:: Attributes Summary

    .. autosummary::
        ~EulerAngles.alpha
        ~EulerAngles.beta
        ~EulerAngles.gamma

    .. rubric:: Methods Summary

    .. autosummary::
        ~EulerAngles.dumps

    .. rubric:: Attributes Documentation
    .. autoattribute:: alpha
    .. autoattribute:: beta
    .. autoattribute:: gamma

    .. rubric:: Method Documentation
    .. automethod:: dumps
