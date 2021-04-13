
.. _understanding_system:

.. .. image:: https://mybinder.org/badge_logo.svg
..  :target: https://mybinder.org/v2/gh/DeepanshS/mrsimulator/master?filepath=jupyternotebooks%2F

.. =====================
.. Optimizing Simulation
.. =====================

.. As developer of ``Mrsimulator`` we put in considerable effort in optimizing the library
.. so that our users have a.

.. Beside the optimization of the underlying C library, the
.. performance of the simulation also depends on the how the user sets the problem. In this
.. section, we will discuss a few points that might speed up simulation.

Understanding isotopomer and method
-----------------------------------
In the previous two sections, we illustrate how to set up a :class:`~mrsimulator.Simulator`
object and run simulations as a three-step process---add spin systems,
add methods, and run the simulation. The Simulator object is designed to handle
bulk simulations and hides most steps from the end-user. One of those steps is how
a method object communications with an isotopomer object to produce a simulation?
The short answer is---each method communicates differently with different spin systems.
You may even consider a method and an isotopomer pair as a complete problem in itself,
and the Simulator object solves multiple such problems.

In this section, we focus our discussion on how a method object communicates with the
spin systems. To understand this, you may need to refresh your NMR 101, may I suggest
chapter 2 of `Understanding NMR spectroscopy` by James Keeler.

Consider the following Bloch decay spectrum method,

.. doctest::

    >>> from mrsimulator.methods import BlochDecaySpectrum
    >>> Bloch_decay_method = BlochDecaySpectrum(
    ...     channels=["13C"], spectral_dimensions=[{"count": 1024, "spectral_width": 25000}]
    ... )

and the following spin systems. A single site spin-1/2 isotopomer is a good
place to start.

**Single spin-1/2 site**

From the introductory quantum mechanic, we know that a single spin-1/2 system is a two
energy levels spin system, where the two energy states are represented as :math:`|\alpha⟩`
and :math:`|\beta⟩` following the bra-ket notation. Here, the energy states
:math:`|\alpha⟩` and :math:`|\beta⟩` corresponds to the quantum numbers :math:`m=+1/2`
and :math:`m=-1/2`, respectively.

Let's create an isotopomer with a single spin-1/2, :math:`^{13}\text{C}` site.

.. doctest::

    >>> from mrsimulator import SpinSystem, Site
    >>> C13 = Site(isotope="13C")
    >>> one_site = SpinSystem(sites=[C13])

To list the Zeeman energy states of an isotopomer, use the
:meth:`~mrsimulator.SpinSystem.zeeman_energy_states` attribute of the instance, as

.. doctest::

    >>> one_site.zeeman_energy_states()
    [|-0.5⟩, |0.5⟩]

A spin system with :math:`n` energy states has :math:`^nP_2` total transitions
connecting any two states, and :math:`n` population. In a single spin-1/2 case,
this accounts to 2 transitions and 2 populations. The two transitions, in bra-ket
notation, follows

.. math::
    |\alpha⟩ \rightarrow |\beta⟩ => |\beta⟩ ⟨\alpha|
    |\beta⟩ \rightarrow |\alpha⟩ => |\alpha⟩ ⟨\beta|

A given method queries the isotopomer for all transitions that are relevant to the
given method. In the case of a BlochDecaySpectrum, the method queries for all
transition with :math:`p = \Delta m = -1`. In the above example, this corresponds to
the :math:`|\beta⟩ ⟨\alpha|` transition. You may access the list of transition
pathways through the :meth:`~mrsimulator.Method.get_transition_pathways` function of
the method object, as

.. doctest::

    >>> relevant_pathways = Bloch_decay_method.get_transition_pathways(one_site)
    >>> print(relevant_pathways)
    [|-0.5⟩⟨0.5|]


**Two spin-1/2 sites**

A two spin-1/2 system is a four energy level spin system, where the energy state are
represented as :math:`|\alpha \alpha⟩`, :math:`|\alpha \beta⟩`, :math:`|\beta \alpha⟩`,
and :math:`|\beta \beta⟩`. There are a total of 12 transitions and 4 populations.
