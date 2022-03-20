.. _query_doc:

Query Objects
=============

There are two types of query objects in mrsimulator---TransitionQuery and MixingQuery.
A TransitionQuery object contains attributes that query the spin systems for
allowed transitions, `i.e.`, transitions that satisfy the query criterion.
A MixingQuery, similarly, queries consecutive transitions for allowed
transition pathways---pathways that satisfy the mixing query criterion.

TransitionQuery
---------------

Transition queries are assigned to the `transition_query` attribute of a SpectralEvent
object as

.. code-block:: python

    from mrsimulator.method import SpectralEvent

    SpectralEvent(
        # other attributes
        transition_query=[
            # list of TransitionQuery objects
        ]
    )

A generic syntax for the :class:`~mrsimulator.method.query.TransitionQuery` object in
python dictionary representation follows,

.. code-block:: python

    t_query = {
        "ch1": {"P": [-1, 1], "D": [0, 2]},
        "ch2": {"P": [+1], "D": [0, -2]},
        "ch3": {"P": [0], "D": [0]},
    }

where ``ch-`` s are the channels over which the query is performed, and its value is a
python dictionary of the :class:`~mrsimulator.method.query.SymmetryQuery` object.
SymmetryQuery is based on the :ref:`spin transition functions <spin_transition_theory>`.
Its attributes,

.. math::
    :label: P_and_D

    P = [(m_f - m_i)_j]_{j=0}^{N_\text{ch}}, \\
    D = [(m_f^2 - m_i^2)_j]_{j=0}^{N_\text{ch}},

are the site-wise list of :math:`p = m_f - m_i` and :math:`d = m_f^2 - m_i^2` transition
symmetry functions, respectively. In Eq. :eq:`P_and_D`, the index :math:`j` runs over the
sites within channel ``ch-`` and :math:`N_\text{ch}` is the total number of sites in the
respective channel. The symbols :math:`m_f` and :math:`m_i` are the spin quantum numbers
for the final and initial energy states of the transition associated with the site at
index :math:`j`. For example, consider the following heteronuclear two-site double-quantum
transition for a two-site spin system [1H, 13C], with ch1=1H and ch2=13C

.. math::
    |0.5, 0.5\rangle \rightarrow |-0.5, -0.5\rangle.

To select this transition, we set up the following query,

.. code-block:: python

    t_query = {
        "ch1": {"P": [-1]},
        "ch2": {"P": [-1]},
    }

The above query selects all transitions where precisely one site in ch1 and one site in
ch2 simultaneously undergoes a :math:`p=-1` transition. Similarly, for spin system [1H, 1H]
with ch1=1H, undergoing a homonuclear two-site double-quantum transition, we write
the query as,

.. code-block:: python

    t_query = {
        "ch1": {"P": [-1, -1]},
    }

Here, the query selects all transitions where exactly two sites in ch1 simultaneously
undergo a :math:`p=-1` transition.

A transition query selects a set of transitions. In the above two examples, the set
consists of a single transition. In the case of a three-site spin system [1H, 1H, 1H],
the same query selects six transitions, where two sites in ch1 simultaneously undergo a
:math:`p=-1` transition, while the third site is at :math:`p=0`.

Note, while the third site defaults to a :math:`p=0`, we do not specify the query as
``"P": [-1, -1, 0]``. The query ``"P": [-1, -1]`` is not the same as ``"P": [-1, -1, 0]``.
While both queries select homonuclear two-site double-quantum transitions, the latter
requires a minimum of a three-site spin system, contrary to a minimum of a two-site spin
system for the first.

Another example of a query object for selecting homonuclear two-site zero-quantum
transitions is

.. code-block:: python

    t_query = {
        "ch1": {"P": [-1, +1]},
    }

In contrast, a query for selecting single-site multi-quantum transitions follow,

.. code-block:: python

    t_query = {
        "ch1": {"P": [-3]},
    }

In the case of a single-site spin system [27Al], the above query will select three
triple-quantum transitions,

.. math::
    |2.5\rangle \rightarrow |-0.5\rangle, \\
    |1.5\rangle \rightarrow |-1.5\rangle, \\
    |0.5\rangle \rightarrow |-2.5\rangle, \\

one symmetric, and two asymmetric transitions. To select the symmetric transition,
modify the query to

.. code-block:: python

    t_query = {
        "ch1": {"P": [-3], "D": [0]},
    }

which first selects the three :math:`p=-3` transitions, and then filters the selection
to transitions where :math:`m_f^2 - m_i^2=0`, `i.e.`,
:math:`|1.5\rangle \rightarrow |-1.5\rangle` central transition.

Rule of Union and Intersection
''''''''''''''''''''''''''''''

As a general rule, the more query criteria we add to the query objects, the smaller
the set of selected transitions. For example, the query from ``"ch1": {"P": [-3]}`` to
``"ch1": {"P": [-3], "D": [0]}`` narrows the transitions from three to one in the case
of a single-site spin system [27Al]. It follows the **intersection** rule---a set
common to all selection criteria.

Now consider a case where we want to select both :math:`p=-1`` and :math:`p=+1`` transitions
simultaneously. Following the rule of intersection, there are precisely zero transitions that
are both :math:`p=+1` and :math:`p=-1`. Here, we use the **union** rule. Recall that the
value of the `transition_query` attribute of the SpectralEvent object is a list of queries,

.. code-block:: python

    SpectralEvent(
        # other attributes
        transition_query=[
            # TransitionQuery(...),  # 0
            # TransitionQuery(...),  # 1
            # TransitionQuery(...),  # 2
        ]
    )

The union rule applies to a set of transitions from multiple transition queries.  In the above
example, the resulting set of selected transitions is the union of transition sets from the
three queries. To select :math:`p=\pm1` transitions, we write

.. code-block:: python

    SpectralEvent(
        # other attributes
        transition_query=[
            # union of set of transitions from query-1 and query-2
            {"ch1": {"P": [-1]}},  # query-1
            {"ch1": {"P": [+1]}},  # query-2
        ]
    )
