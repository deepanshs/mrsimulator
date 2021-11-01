.. _installing_python

Installing Python
-----------------

.. warning::
  Some of the dependency libraries for ``mrsimulator`` are incompatible with Python 3.10.0.
  We anticipate the issues with Python 3.10.0 to be fixed in the coming months.

  Until then, we recommend installing ``mrsimulator`` in a Python 3.9 environment. See the
  section on setting up virtual environments.

Installing Python
"""""""""""""""""

``mrsimulator`` requires Python or a hosted Notebook service to run. If you are
using Google Colab, see the `colab instillation steps <on_google_colab>`__ or the
`colab issues <Google Colab Notebook Issues_>`__.

To check if python is installed on your system, open a new terminal/command line and run
**Checking the version of Python**

.. tabs::

    .. tab:: Windows

        To check if Python is installed on a Windows machine, first open the Command Prompt application.
        Next, type

        .. code-block:: shell

            python -V

        and press enter. ``Python 2.7.x`` or ``Python 3.x.x`` should be printed to the console.
        However, if Python is not installed, an error message like the following will appear:

        .. code-block:: shell

            'python' is not recognized as an internal or external command, operable program or batch file.

        To install Python, visit `python.org <https://www.python.org/downloads/>`__ to download a
        version of Python 3. During the instillation process, please check **add Python 3.x to PATH** box.

        If you are certain Python is installed on your system but continue to receive errors,
        `check if Python is in your PATH variable <https://datatofish.com/add-python-to-windows-path/>`__.

    .. tab:: MacOS

        Note that most recent versions of MacOS come with Python pre-installed. If you're unsure
        if Python is installed, follow these steps.

        To check if Python is installed on MacOS, open the Terminal application. Next, type

        .. code-block:: bash

            python -V

        and press enter. ``Python 2.7.x`` or ``Python 3.x.x`` should be printed to the console.
        However, if Python is not installed, an error message like the following will appear:

        .. code-block:: shell

            zsh: command not found: python

        To install Python, visit `python.org <https://www.python.org/downloads/>`__ to download a
        version of Python 3 for your system.

    .. tab:: Linux

        Note that most recent versions of MacOS come with Python pre-installed. If you're unsure
        if Python is installed, follow these steps.

        To check if Python is installed on Linux, open a terminal. Next, type

        .. code-block:: bash

            python -V

        and press enter. ``Python 2.7.x`` or ``Python 3.x.x`` should be printed to the console.
        However, if Python is not installed, an error message like the following will appear:

        .. code-block:: shell

            bash: python: command not found

        To install Python, visit `python.org <https://www.python.org/downloads/>`__ to download a
        version of Python 3 for your system.


Updating Python
"""""""""""""""

If Python is already installed on your system but is out of date, we recommend `installing Anaconda
<Virtual Environments for Python>`__ to manage Python versions. Anaconda is versatile and allows
multiple versions of Python to run on one computer without interfering with each-other.

However, if Anaconda can't be used, newer versions of Python can be installed from `python.org
<https://www.python.org/downloads/>`__. We recommend using the latest version of ``Python 3.9``.
