
PIP is a package manager for Python packages and is included with python version 3.4
and higher. PIP is the easiest way to install python packages. Install the package
using pip as follows,

.. code-block:: bash

    $ pip install mrsimulator

For *Mac* users, if the above statement didn't work, you are probably using mac OS
system python, in which case, use the following,

.. code-block:: bash

    $ python3 -m pip install mrsimulator --user

If you get a ``PermissionError``, it usually means that you do not have the required
administrative access to install new packages to your Python installation. In this
case, you may consider adding the ``--user`` option at the end of the statement to
install the package into your home directory. You can read more about how to do this in
the `pip documentation <https://pip.pypa.io/en/stable/user_guide/#user-installs>`_.
