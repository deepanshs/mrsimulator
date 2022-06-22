.. _install_python:

Installing Python
"""""""""""""""""

``mrsimulator`` requires Python or a hosted Notebook service to run. If you are
using Google Colab, see specific instructions in the :ref:`install` section.

Checking the version of Python
''''''''''''''''''''''''''''''

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
        version of Python 3. During the installation process, check **add Python 3.x to PATH**.
        If this isn't selected, Python may not be accessible across your computer and will cause
        errors.

        If you are certain Python is installed on your system but continue to receive errors, see
        `adding Python to your PATH variable <https://datatofish.com/add-python-to-windows-path/>`__.

    .. tab:: MacOS

        Most recent versions of MacOS come with Python pre-installed. If you're unsure
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
'''''''''''''''

If Python is already installed on your system but is out of date, we recommend `installing Anaconda
<Virtual Environments for Python>`__ to manage Python versions. Anaconda is versatile and allows
multiple versions of Python to run on one computer without interfering with each-other.

However, if Anaconda can't be used, newer versions of Python can be installed from `python.org
<https://www.python.org/downloads/>`__. We recommend using the latest version of ``Python 3.9``.
