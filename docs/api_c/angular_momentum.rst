


Angular Momentum Method Documentation
-------------------------------------

Generic methods
^^^^^^^^^^^^^^^

.. doxygenfunction:: wigner_d_element()

.. doxygenfunction:: wigner_d_matrices()


.. Specialized methods
.. ^^^^^^^^^^^^^^^^^^^
.. Evaluating sines and cosines of an angle is computationally expensive.
.. When the cosine and sine of the angles are known, we provide specialized
.. angular momentum methods where the argument is :math:`\exp(\theta)`, expressed
.. as :math:`\cos(\theta) + I \sin(\theta)`.

.. Within the methods, the :math:`\exp(\theta)` is stored as an array of length
.. two, where the first element is :math:`\cos(\theta)` and, the second element is
.. :math:`\sin(\theta)`.

.. .. doxygenfunction:: wigner_d_element_from_exp_I_beta()

.. .. doxygenfunction:: wigner_d_matrices_from_exp_I_beta()

.. .. doxygenfunction:: __wigner_rotation_2()
