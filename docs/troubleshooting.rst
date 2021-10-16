.. _troubleshooting:

###########################
Common Issues and Solutions
###########################

.. _instillation_issues:

Instillation Issues
-------------------

**Installing Python**

``mrsimulator`` requires Python or a hosted Notebook service to run. If you are
using Google Colab, see the `colab instillation steps <on_google_colab>`__ or the
`colab issues <Google Colab Notebook Issues_>`__.

To check if python is installed on your system, open a new terminal/command line and run

.. code-block:: shell

      python -V

which should output ``Python 2.x.x`` or ``Python 3.x.x`` if Python is installed. Otherwise, visit
`python.org <https://www.python.org/downloads/>`__ to install Python.

.. warning::
  ``mrsimulator`` is incompatible with Python 3.10.0 (released October 4th, 2021) due to issues 
  when installing dependencies. We anticipate the issues with Python 3.10.0 to be fixed in 
  the coming months.
  
  If ``python -V`` prints ``3.10.0``, please follow instructions in
  :ref:`virtual_envs` to install an older version of python.

**Installing PIP**

PIP is the dominant python package manager. If you encounter a message like

.. code-block:: shell

      'pip' is not recognized as an internal or external command, operable program or batch file.

please follow the `pip instillation instructions <https://pip.pypa.io/en/stable/installation/>`__
to install pip on your system. Afterwards continue following the :ref:`instillation instructions <install>`.


Google Colab Notebook Issues
----------------------------

Google Colab is an extremely accessible platform for running Python code, but some issues can arise
when using ``mrsimulator`` on google colab

**Updating numpy in Google Colab**

If you see the error

.. code-block:: bash

      ValueError: numpy.ndarray size changed, may indicate binary incompatibility. Expected 88 from C header, got 80 from PyObject

that means an older version of ``numpy`` is installed. To update numpy to a compatable version,
type 

.. code-block:: shell

      !pip install -U numpy

in a new cell and run. The runtime for Google Colab needs to be restarted to use the newer version.
After restarting the runtime, everything ``numpy`` should run without error.


.. _virtual_envs:

Virtual Environments for Python
-------------------------------

**Creating a Python environment using Anaconda**

Since different Python libraries have different dependencies, installing multiple libraries in the
same place can cause issues. For example, ``mrsimulator`` requires at least ``numpy v1.17`` but
``some-other-library`` might require exactly ``numpy v1.15``. These two libraries would likely
throw errors when run in the same environment.

For this reason, we recommend using an environment manager, like ``vnev`` or ``anaconda``.
We will look at anaconda for its simple commands. Instillation instructions can be found on the
`anaconda documentation page <https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html>`__.

.. note:: 
  Anaconda is a robust package and environment management program, but it does require
  a significant amount of space on disk (>400mb). If you need a lightweight environment manager
  and are confident with python we recommend looking at `Python's venv documentation
  <https://docs.python.org/3/library/venv.html>`__.

Once Anaconda is installed, create a new environment by running

.. code-block:: shell

      conda create -n <name> python=3.9

where ``<name>`` is the desired name of your environment. Next activate the environment by running

.. code-block:: shell
      
      conda activate <name>

To install ``mrsimulator`` in this new environment, follow the :ref:`instillation <install>` 
instructions. To exit the environment run

.. code-block:: shell

      conda deactivate

To start using ``mrsimulator`` again, simply activate the same environment.

Packages installed in an environment remain installed between sessions and won't interfere
with packages in other environments.


Still Having Issues?
--------------------

If this page didn't resolve your issue or if you are having problems with ``mrsimulator``
please submit an bug report on our `Github issue tracker <https://github.com/DeepanshS/mrsimulator/issues>`_.

Github is the preferred bug and issue reporting method, but issues can also be reported by
directly contacting `Matthew Giammar <mailto:giammar.7@osu.edu>`__. 

Discussions are welcome on our `Github discussion <https://github.com/DeepanshS/mrsimulator/discussions>`_
page.
