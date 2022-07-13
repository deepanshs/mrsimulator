.. _mixing_query_documentation:

==============
Mixing Queries
==============

The amplitude of a transition pathway signal derives from the product
of mixing amplitudes associated with each transfer between transitions in a
transition pathway,
e.g.,

.. math::
    (a_{0A})\,\hat{t}_A \rightarrow (a_{0A}a_{AB})\,\hat{t}_B\rightarrow (a_{0A}a_{AB}a_{BC})\,\hat{t}_C \rightarrow \cdots

Here, :math:`a_{0A}` is the amplitude of the initial :math:`\hat{t}_A` transition,
:math:`a_{AB}` is the mixing amplitude for the transfer from
:math:`\hat{t}_A` to :math:`\hat{t}_B`,  :math:`a_{BC}` is the mixing amplitude
for the transfer from :math:`\hat{t}_B` to :math:`\hat{t}_C`, and so on.  The
growing product :math:`(a_{0A}a_{AB}a_{BC} \cdots)` is the transition pathway
amplitude.  Eliminating a transition with a TransitionQuery in a spectral or
delay event sets the eliminated transition's pathway amplitude to zero, i.e., it
prunes that transition pathway branch.

Default Total Mixing between Adjacent Spectral or Delay Events
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

In previous discussions, we did not mention the efficiency of transfer between
selected transitions in adjacent **SpectralEvent** objects. This is because, as
default behavior, **mrsimulator** does a *total mixing*, i.e., connects all
selected transitions in the two adjacent spectral or delay events. In other
words, if the first of two adjacent **SpectralEvent** objects has three selected
transitions, and the second has two selected transitions, then **mrsimulator**
will make :math:`3 \times 2 = 6` connections, i.e., six transition pathways
passing from the first to second **SpectralEvent** objects.

Additionally, this total mixing event assumes that every connection has a mixing
amplitude of 1. This is unrealistic, but if used correctly gives a significant
speed-up in the simulation by avoiding the need to calculate mixing amplitudes.

.. warning::

    The use of total mixing, i.e., the default mixing, can complicate the comparison
    of integrated intensities between different methods, depending on the selected
    transition pathways.

If this default mixing behavior had been explicitly shown in the previous example,
the events list in the first **SpectralDimension** would have looked like the code
below.

.. plot::
    :context: close-figs

    from mrsimulator import Site, SpinSystem, Simulator
    from mrsimulator.method import SpectralEvent, MixingEvent

    events = [
        SpectralEvent(fraction=9 / 16, transition_queries=[{"ch1": {"P": [-3], "D": [0]}}]),
        MixingEvent(query="TotalMixing"),
        SpectralEvent(fraction=7 / 16, transition_queries=[{"ch1": {"P": [-1], "D": [0]}}]),
        MixingEvent(query="TotalMixing"),
    ]

Since only one transition was selected in each **SpectralEvent**, the expected (and
default) behavior is that there is a mixing (transfer) of coherence between
the symmetric triple-quantum and central transitions, forming the desired
transition pathway.

However, when multiple transition pathways are present in a method, you may need
more accurate mixing amplitudes when connecting selected transitions of adjacent
events. You may also need to prevent the undesired mixing of specific
transitions between two adjacent events. As described below, you can avoid a
``"TotalMixing"`` event by inserting **MixingEvent** object with a certain rotation
query.

Rotation Query
''''''''''''''

A rotation of :math:`\theta` about an axis defined by :math:`\phi`  in the
:math:`x`-:math:`y` plane on a selected transition, :math:`\ketbra{I, m_f}{I,
m_i}`, in a spectral or delay event transfers it to all selected transitions,
:math:`\ketbra{I,m_f'}{I,m_i'}` in the next spectral or delay event, according
to

.. math::

    \ketbra{I, m_f}{I, m_i} \stackrel{\theta_\phi}{\longrightarrow} \sum_{m_f'}\sum_{m_i'} d_{m_f',m_f}^{(I)}(\theta)d_{m_i',m_i}^{(I)}(\theta)e^{-i\Delta p\phi}(i)^{\Delta p}\ketbra{I,m_f'}{I,m_i'},

where :math:`\Delta p_I = p_I' - p_I`.  From this expression, we obtain the complex mixing amplitude
from :math:`\ketbra{I, m_f}{I, m_i}` to :math:`\ketbra{I, m_f'}{I, m_i'}` due to a rotation to be

.. math::

    a(\theta,\phi) = d_{m_f',m_f}^{(I)}(\theta)d_{m_i',m_i}^{(I)}(\theta)e^{-i\Delta p\phi}(i)^{\Delta p}.

From this expression, we note a few interesting and useful cases.  One is the coherence
transfer under a :math:`\pi` rotation, given by

.. math::
    :label: piPulseTransition

    \ketbra{I,m_f}{I, m_i}  \stackrel{\pi_\phi}{\longrightarrow} \ketbra{I, -m_f}{I, -m_i} e^{-i\Delta p\phi}(i)^{\Delta p}.

That is, a :math:`\pi` rotation will make only one connection between
transitions in adjacent events.  It is also a special connection because the
:math:`\text{p}_I` transition symmetry value for the two transitions are equal
but opposite in sign.  Additionally, the :math:`\text{d}_I` transition symmetry
remains unchanged (:math:`\Delta \text{d}_I = 0`) for the two transitions. (In
fact, this behavior under a :math:`\pi` rotation is generally true for odd
(:math:`\text{p}_I, \text{f}_I, \ldots)` and even (:math:`\text{d}_I,
\text{g}_I, \ldots)` rank spin transition symmetry functions.)

Another interesting result is that, while a rotation can transfer a transition
into many other transitions, the :math:`\text{d}_I` transition symmetry value
cannot remain unchanged (:math:`\Delta \text{d}_I \neq 0`) between two connected
transitions under a :math:`\pi/2` rotation.

Finally, another useful result is

.. math::
    :label: zeroPulseTransition

    \ketbra{I,m_f}{I, m_i}  \stackrel{(0)_\phi}{\longrightarrow} \ketbra{I, m_f}{I, m_i}.

While it's not surprising that a rotation through an angle of zero does nothing
to the transition, this turns out to help act as the opposite of a total mixing
event, i.e., a **no mixing** event. As a convenience, this is defined as a
``"NoMixing"`` query and can be implemented with the code below.

.. plot::
    :context: close-figs

    MixingEvent(query="NoMixing")

The **MixingEvent** object holds the rotation details in a **MixingQuery** object as
a **RotationQuery** object associated with a ``channels`` attribute.  This is
illustrated in the sample code below.

.. plot::
    :context: close-figs

    import numpy as np
    from mrsimulator.method.query import RotationQuery
    rot_query_90 = RotationQuery(angle=np.pi/2, phase=0)
    rot_query_180 = RotationQuery(angle=np.pi, phase=0)
    rot_mixing = MixingEvent(query={
            "ch1": rot_query_90,
            "ch2": rot_query_180
        }
    )
