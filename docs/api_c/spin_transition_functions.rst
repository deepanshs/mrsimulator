.. _spin_transition_function:

Spin transition functions (STF), :math:`\xi_L^{(k)}(i,j)`
---------------------------------------------------------

.. seealso:: :ref:`Spin transition functions <spin_transition_theory>`

.. raw:: html

    <a class="btn btn-default"
       href=./source/spin_transition_functions.html> Source
    </a>


Single nucleus spin transition functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The single spin transition functions for
:math:`\left|m_i\right> \rightarrow \left|m_f\right>` transition, where :math:`m_j`
is the spin quantum number, and the subscripts :math:`i` and :math:`f` refer to the
initial and final energy states.

.. doxygenfunction:: STF_p
   :project: mrsim

.. doxygenfunction:: STF_d(const double, const double)
   :project: mrsim

.. doxygenfunction:: STF_f
   :project: mrsim

Composite single nucleus spin transition functions
""""""""""""""""""""""""""""""""""""""""""""""""""

The composite single spin transition functions are linear combinations of the
single spin transition functions.

.. doxygenfunction:: STF_cL
   :project: mrsim

Two weakly coupled nuclei spin transition functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The weakly coupled spin transition function for
:math:`\left|m_{i_I}, m_{i_S}\right> \rightarrow \left|m_{f_I}, m_{f_S}\right>`
transition. Here, the subscript :math:`I` and :math:`S` denotes the two weakly
coupled spins.

.. doxygenfunction:: STF_dIS
   :project: mrsim
