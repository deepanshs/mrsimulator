
Colaboratory is a Google research project. It is a Jupyter notebook environment that
runs entirely in the cloud. Launch a new notebook on
`Colab <https://colab.research.google.com>`_. We recommend going through the *Welcome to Colab!*
tutorial if you are new to Notebooks.

.. note::
    Check the version of ``numpy`` on colab by typing ``!pip show numpy`` in the first cell. If
    the version less that 1.20, type ``!pip install numpy -U`` in the next cell to update numpy.
    Once the update finishes, restart the kernel by pressing *Runtime -> Restart Runtime* button.

To install the ``mrsimulator`` package, type

.. code-block:: shell

    !pip install mrsimulator

in a new cell, and execute. All done! You may now start using the library, or
proceed to :ref:`getting_started` to continue the tutorial.
