.. _troubleshooting:

###########################
Common Issues and Solutions
###########################

.. _installation_issues:

Installation Issues
-------------------

**Installing Python**

``mrsimulator`` requires Python or a hosted Notebook service to run. If you are
using Google Colab, see the `colab installation steps <on_google_colab>`__ or the
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

please follow the `pip installation instructions <https://pip.pypa.io/en/stable/installation/>`__
to install pip on your system. Afterwards continue following the :ref:`installation instructions <install>`.



.. _virtual_envs:




Still Having Issues?
--------------------

If this page didn't resolve your issue or if you are having problems with ``mrsimulator``
please submit an bug report on our `Github issue tracker <https://github.com/DeepanshS/mrsimulator/issues>`_.

Github is the preferred bug and issue reporting method, but issues can also be reported by
directly contacting `Matthew Giammar <mailto:giammar.7@osu.edu>`__.

Discussions are welcome on our `Github discussion <https://github.com/DeepanshS/mrsimulator/discussions>`_
page.
