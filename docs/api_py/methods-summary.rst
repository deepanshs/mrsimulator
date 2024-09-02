.. _methods_summary_api:

Methods
=======

.. currentmodule:: mrsimulator.method.lib

The following are the list of methods currently supported by **MRSimulator** as a part
of the ``mrsimulator.method.lib`` module. To import a method, for example the
*BlochDecaySpectrum*, used

.. doctest::

    >>> from mrsimulator.method.lib import BlochDecaySpectrum

.. All methods categorize into two groups, generic and specialized methods. A generic
.. method is general and is based on the number of spectral dimensions. At present,
.. there are two generic methods, ``Method1D`` and ``Method2D``. All specialized methods
.. are derived from their respective generic method objects. The purpose of the specialized
.. methods is to facilitate user ease when setting up some commonly used methods, such as
.. the MQMAS, STMAS, PASS, MAT, etc.


Summary
-------

.. **Generic methods**
..
.. .. autosummary::
..     ~Method1D
..     ~Method2D

**Specialized 1D methods**

.. autosummary::
    ~BlochDecaySpectrum
    ~BlochDecayCTSpectrum

**Specialized 2D methods**

.. autosummary::
    ~ThreeQ_VAS
    ~FiveQ_VAS
    ~SevenQ_VAS
    ~ST1_VAS
    ~ST2_VAS
    ~SSB2D

.. ~Cosy

.. UML Diagram
.. -----------

.. .. figure:: ../_static/classes_methods.*


Table of contents
-----------------

.. toctree::
    methods/BlochDecaySpectrum
    methods/BlochDecayCTSpectrum
    methods/mqvas_api
    methods/stvas_api
    methods/ssb2d_api

.. methods/cosy
