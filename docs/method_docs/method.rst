
.. _method_documentation:

======
Method
======

While **mrsimulator**'s organization of the :ref:`spin_sys_api` object and its
composite objects, :ref:`site_api`, and :ref:`coupling_api` are easily
understood by anyone familiar with the underlying physical concepts, the
organization of the :ref:`method_api` object in **mrsimulator** and its related
composite objects require a more detailed explanation of their design. This
section assumes that you are already familiar with the topics covered in the
Introduction sections :ref:`getting_started`,
:ref:`introduction_isotopomers_example`, and :ref:`fitting_example`.

.. note::

    Lists in Python are always ordered; however, we use the phrase *ordered list* when
    referring to a list where the order of elements is meaningful for **mrsimulator**.
    Similarly, the phrase *unordered list* implies **mrsimulator** does not care about
    the ordering of the items.

Overview
--------

An experimental NMR method involves a sequence of rf pulses, free evolution
periods, and sample motion. The **Method** object in **mrsimulator** models the
spectrum from an NMR pulse sequence. The **Method** object is designed to be
versatile in its ability to model spectra from various multi-pulse NMR methods
using concepts from the `symmetry pathway approach
<https://doi.org/10.1016/j.pnmrs.2010.11.003>`_ where a pulse sequence is
understood in terms of a set of desired (and undesired)  *transition pathways*.
Each transition pathway is associated with a single resonance in a
multi-dimensional NMR spectrum. The transition pathway signal encodes
information about the spin system interactions in its amplitude and correlated
frequencies. Consider the illustration of a 2D pulse sequence shown below, where
a desired signal for the method is associated with a particular transition
pathway, :math:`{\hat{A} \rightarrow \hat{B} \rightarrow \hat{C} \rightarrow
\hat{D} \rightarrow \hat {E} \rightarrow \hat{F}}`.

.. figure:: ../_static/TransitionPathway.*
    :width: 600
    :alt: figure
    :align: center

Here, the first spectral dimension, i.e., the Fourier transform of the
transition pathway signal as a function of :math:`t_1`, derives its *average
frequency*, :math:`\overline{\Omega}_1`, from a weighted average of the
:math:`\hat{A}`, :math:`\hat{B}`, and :math:`\hat{C}` transition frequencies.
The second spectral dimension, i.e., the FT with respect to :math:`t_2`, derives
its average frequency, :math:`\overline{\Omega}_2`, from a weighted average of
the :math:`\hat{E}`, and :math:`\hat{F}` transition frequencies. Much of the
experimental design and implementation of an NMR method is in identifying the
desired transition pathways and finding ways to acquire their signals while
eliminating all undesired transition pathway signals.

While NMR measurements occur in the time domain, **mrsimulator** simulates the
corresponding multi-dimensional spectra directly in the frequency domain. The
**Method** object in **mrsimulator** needs only a few details of the NMR pulse
sequence to generate the spectrum. It mimics the result of the pulse sequence
given the desired transition pathways and their complex amplitudes and average
frequencies in each spectroscopic dimension of the dataset. To this end, a
**Method** object is organized according to the UML diagram below.

.. figure:: ../_static/MethodUML.*
    :width: 700
    :alt: figure
    :align: center

.. note::

 In UML (Unified Modeling Language) diagrams, each class is represented with a
 box that contains two compartments. The top compartment has the class's name,
 and the bottom compartment contains the class's attributes. Default attribute
 values are shown as assignments. A composition is depicted as a binary
 association decorated with a filled black diamond. Inheritance is shown as a
 line with a hollow triangle as an arrowhead.

At the heart of a **Method** object, assigned to its attribute
``spectral_dimensions``, is an ordered list of :ref:`spectral_dim_api` objects
in the same order as the time evolution dimensions of the experimental NMR
sequence. In each **SpectralDimension** object, assigned to the attribute
``events``, is an ordered list of :ref:`event_api` objects, which are divided
into three types: (1) :py:meth:`~mrsimulator.method.SpectralEvent`, (2)
:py:meth:`~mrsimulator.method.DelayEvent`, and (3)
:py:meth:`~mrsimulator.method.MixingEvent`.  This ordered list of Event objects
is used to select the desired transition pathways and determine their average
frequency and complex amplitude in the **SpectralDimension**.

.. warning::

  DelayEvent objects are not available in version 0.7 of **mrsimulator**.

**SpectralEvent** and **DelayEvent** objects define which transitions are
observed during the event and under which transition-dependent frequency
contributions they evolve. No coherence transfer among transitions or
populations occurs in a spectral or delay event. The transition-dependent
frequency contributions during an Event are selected from a list of
:ref:`enumeration literals<freq_contrib_api>` and placed in the ``freq_contrib``
attribute of the event. If ``freq_contrib`` is left unspecified, i.e., the
value of ``freq_contrib`` is set to ``None``, a default list holding the
enumeration literals for *all* contributions is generated for the event.

.. note::

  All frequency contributions from direct and indirect spin-spin couplings are
  calculated in the weak-coupling limit in **mrsimulator**.

Additionally, the user can affect transition frequencies during a spectral or
delay event by changing other measurement attributes: ``rotor_frequency``,
``rotor_angle``, and ``magnetic_flux_density``. If left unspecified, these
attributes default to the values of the identically named global attributes in
the **Method** object. **SpectralEvent** objects use the ``fraction`` attribute to
calculate the weighted average frequency during the spectral dimension for each
selected transition pathway.

Inside **SpectralEvent** and **DelayEvent** objects, is a list of
:py:meth:`~mrsimulator.method.query.TransitionQuery` objects (*vide infra*)
which determine which transitions are observed during the event. **Method**
objects in **mrsimulator** are general-purpose because they are designed for an
arbitrary spin system, i.e., a method does not know the spin system in advance.
When designing a **Method** object, you cannot identify and select a transition
through its initial and final eigenstate quantum numbers. Transition selection
is done through **TransitionQuery** and
:py:meth:`~mrsimulator.method.query.SymmetryQuery` objects during individual
spectral or delay events. **TransitionQuery** objects can hold a
**SymmetryQuery** object in ``ch1``, ``ch2``, or ``ch3``, which act on
specific isotopes defined by the ``channels`` attribute in **Method**. It is
only during a simulation that the **Method** object uses its **TransitionQuery**
objects to determine the selected transition pathways for a given **SpinSystem**
object by the initial and final eigenstate quantum numbers of each transition.

Between adjacent SpectralEvent or DelayEvent objects, **mrsimulator** defaults
to *total mixing*, i.e., connecting all selected transitions in the two adjacent
spectral or delay events. This default behavior can be overridden by placing an
explicit **MixingEvent** object between such events. Inside **MixingEvent**
objects is a :py:meth:`~mrsimulator.method.query.MixingQuery` object, which
determines the coherence transfer amplitude between transitions. A
**MixingQuery** object holds
:py:meth:`~mrsimulator.method.query.RotationQuery` objects acting on specific
isotopes in the spin system. As before, the isotope upon which the
**RotationQuery** objects act is determined by the ``channels`` attribute in the
**Method** object.

In this guide to designing custom Method objects, we begin with a brief review
of the relevant *Symmetry Pathway* concepts employed in **mrsimulator**. This
review is necessary for understanding (1) how transitions are selected during
spectral and delay events and (2) how average signal frequencies and amplitudes
in each spectral dimension are determined. We outline the procedures for
designing and creating **TransitionQuery** and **MixingQuery** for single- and
multi-spin transitions and how to use them to select the transition pathways
with the desired frequency and amplitudes in each **SpectralDimension** of your
custom method object. In multi-dimensional spectra, we illustrate how the
desired frequency correlation can sometimes be achieved by using an appropriate
affine transformation. We also examine how changing the frequency contributions
in **SpectralEvent** of **DelayEvent** objects can be used to obtain the desired
frequency and amplitude behavior. The ability to select :ref:`frequency
contributions<freq_contrib_api>` can often reduce the number of events needed in
the design of your custom Method object.


.. Sections
.. --------

.. These sections need to be converted from a toctree to a list of page references

.. .. toctree::
..     :maxdepth: 1

..     transition_query
..     mixing_query
..     frequency_contrib
..     affine_transformation
..     origin_and_reference_offset
..     method_theory
..     attribute_tables
