.. _virtual_environment_troubleshooting

Virtual Environments
--------------------

Creating a Python environment using Anaconda
""""""""""""""""""""""""""""""""""""""""""""

Since different Python libraries have different dependencies, installing multiple libraries in the
same place can cause issues. For example, ``mrsimulator`` requires at least ``numpy v1.17`` but
``some-other-library`` might require exactly ``numpy v1.15``. These two libraries would likely
throw errors when run in the same environment.

For this reason, we recommend using an environment manager, like ``venv`` or ``anaconda``.
We will look at anaconda for its simple commands. Instillation instructions can be found on the
`anaconda documentation page <https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html>`__.

.. note::
  Anaconda is a robust package and environment management program, but it does require
  a significant amount of space on disk (>400mb). If you need a lightweight environment manager
  and are confident with Python, we recommend looking at `Python's venv documentation
  <https://docs.python.org/3/library/venv.html>`__.

Once Anaconda is installed, create a new environment by running

.. code-block:: shell

    $ conda create -n <name> python=3.9

where ``<name>`` is the desired name of your environment. Each environment can have a Python
version specified after ``python=``. We recommend using ``python=3.9``. Next activate the environment by running

.. code-block:: shell

    $ conda activate <name>

The current environment name is reflected in the leftmost portion of a line. On MacOS, an
environment named ``mrsimulator-0.7`` should look like

.. code-block:: shell

    (mrsimulator-0.7) nmruser@machine $

If you are using a code editor or IDE, the current environment should be displayed somewhere on
the window. For VS Code, the environment name and Python version are shown in the bottom-left
corner.

To install ``mrsimulator`` in this new environment, follow the :ref:`instillation <install>`
instructions. ``mrsimulator`` and any other libraries will only be installed in the active
environment. This way different projects can run in separate environments.

To exit the environment run

.. code-block:: shell

    $ conda deactivate

To start using ``mrsimulator`` again, simply activate the environment in which it was installed.

Packages installed in an environment remain installed between sessions and won't interfere
with packages in other environments.
