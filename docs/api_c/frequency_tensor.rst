
.. _frequency_component_function:

Frequency Tensors (FT), :math:`\Lambda_{L, n}^{(k)}(i,j)`
---------------------------------------------------------

.. seealso:: :ref:`Frequency component functions in PAS <frequency_tensor_theory>`,
             :ref:`Spatial orientation functions (SOF) <spatial_orientation_function>`,
             :ref:`Spin transition functions (STF) <spin_transition_function>`

.. raw:: html

    <a class="btn btn-default"
       href=./source/frequency_tensor_source.html> Source
    </a>

Single nucleus frequency tensor components
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First-order Nuclear shielding
"""""""""""""""""""""""""""""

.. doxygenfunction:: FCF_1st_order_nuclear_shielding_tensor_components()

First-order Electric Quadrupole
"""""""""""""""""""""""""""""""

.. doxygenfunction:: FCF_1st_order_electric_quadrupole_tensor_components()

Second-order Electric Quadrupole
""""""""""""""""""""""""""""""""

.. doxygenfunction:: FCF_2nd_order_electric_quadrupole_tensor_components()



Two coupled nucleus frequency tensor components
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First-order J-coupling (weak coupling limit)
""""""""""""""""""""""""""""""""""""""""""""

.. doxygenfunction:: FCF_1st_order_weak_J_coupling_tensor_components()
