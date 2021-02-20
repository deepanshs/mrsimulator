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
method is general and is based on the number of spectral dimensions. At present,
there are two generic methods, ``Method1D`` and ``Method2D``. All specialized methods
are derived from their respective generic method objects. The purpose of the specialized
methods is to facilitate user ease when setting up some commonly used methods, such as
the MQMAS, STMAS, PASS, MAT, etc.


Summary
-------

**Generic methods**

.. autosummary::
    ~Method1D
    ~Method2D

**Specialized methods**

.. autosummary::
    ~BlochDecaySpectrum
    ~BlochDecayCTSpectrum
    ~ThreeQ_VAS
    ~FiveQ_VAS
    ~SevenQ_VAS
    ~ST1_VAS
    ~ST2_VAS
    ~SSB2D


Table of contents
-----------------

.. toctree::
    methods/method1D
    methods/BlochDecaySpectrum
    methods/BlochDecayCTSpectrum
    methods/method2D
    methods/MQVAS
    methods/stvas
    methods/SSB2D
