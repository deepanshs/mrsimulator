For developers and contributors
===============================

Setting up a dedicated code editor
''''''''''''''''''''''''''''''''''

Using a code editor or IDE is useful when contributing to a codebase. There are many products
available and you may use whatever you're comfortable with. For new developers, we recommend
`VS Code <https://code.visualstudio.com>`_ since it is lightweight, free and has a breadth of
community extensions.

Make your own copy of mrsimulator on GitHub
'''''''''''''''''''''''''''''''''''''''''''

Making a copy of someone's code on GitHub is the same as making a *fork*. A fork is a
complete copy of the code and all of its revision history.

1. Log into a `GitHub account <https://github.com>`_.
2. Go to the `mrsimulator Github <https://github.com/deepanshs/mrsimulator>`_ home page.
3. Click on the *fork* button.

You will see a short animation of Octocat scanning a book on a flatbed scanner. After
that, you should find yourself at the home page for your forked copy of mrsimulator.


Create a development environment
''''''''''''''''''''''''''''''''

It is good practice to create separate virtual python environments for packages when
in developing packages. There are many environment managers available, however, we recommend using
`Anaconda or Miniconda <https://docs.anaconda.com/anaconda/install/>`_ for Linux, Mac with Intel
processors, and windows.

.. note::
    *Anaconda Python is **incompatible** with Apple silicon based Macs* Follow the instructions
    under :ref:`OS-dependent prerequisites` for steps to install Arm64 compatible Python.

The following is an example of creating Conda environment

.. code-block:: bash

    $ conda create -n mrsimulator-dev python=3.9

The above command will create a new environment named *mrsimulator-dev* using python 3.9. To
activate the environment, use

.. code-block:: bash

    $ conda activate mrsimulator-dev


Make sure git is installed on your computer
'''''''''''''''''''''''''''''''''''''''''''

`Git <https://git-scm.com>`_ is the name of a source code management system. It keeps
track of the changes made to the code and manages contributions from several different
individuals. You may read more about git at the `Git Basics <https://git-scm.com/book/>`_.

If you are using anaconda/miniconda, you probably have git pre-installed. To check, type
in terminal

.. code-block:: bash

    $ git --version
    # if git is installed, will get something like: git version 2.30.2

If git is not installed, `install <https://git-scm.com/downloads>`_ it.


**Basic git configuration:**

Follow the instructions at `Set Up Git <https://docs.github.com/en/github/getting-started-with-github/set-up-git#set-up-git>`_
at GitHub to configure:

- Your user name and email in your copy of git.
- Authentication, so you don’t have to type your GitHub password every time you need to
  access GitHub from the command line.


Copy your fork of mrsimulator from GitHub to your computer
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Unless you plan on always editing the code using the online Github editor, you may need to
copy the fork of mrsimulator from your GitHub account to your computer. Make a complete
copy of the fork with

.. code-block:: bash

    $ git clone --recursive https://github.com/your-user-name/mrsimulator.git

Insert *your-user-name* with your GitHub account username. If there is an error at this
stage, it is probably an error in setting up authentication.

You now have a copy of the mrsimulator fork from your GitHub account to your local computer
into a mrsimulator folder.

Understanding *Remotes*
'''''''''''''''''''''''

In git, the name for another location of the same repository is *remote*.
The repository that contains the latest "official" development version is traditionally
called the *upstream* remote. You can read more about
`remotes on Git Basics <https://git-scm.com/book/en/v2/Git-Basics-Working-with-Remotes>`_.

At this point, your local copy of mrsimulator doesn't know where the *upstream* development
version of mrsimulator is. To let git know, change into the mrsimulator folder you created in
the previous step, and add a remote:

.. code-block:: bash

    cd mrsimulator
    git remote add mrsimulator git://github.com/deepanshs/mrsimulator.git

You can check that everything is set up properly so far by asking git to show you all of the
remotes it knows about for your local repository of mrsimulator with ``git remote -v``, which
should display something like

.. code-block:: bash

    mrsimulator   git://github.com/deepanshs/mrsimulator.git (fetch)
    mrsimulator   git://github.com/deepanshs/mrsimulator.git (push)
    origin     git@github.com:your-user-name/mrsimulator.git (fetch)
    origin     git@github.com:your-user-name/mrsimulator.git (push)


Build the development version of mrsimulator
''''''''''''''''''''''''''''''''''''''''''''

OS-dependent prerequisites
""""""""""""""""""""""""""

.. note::
    Installing OS-dependent prerequisites is a one-time process. If you are
    upgrading to a newer version of mrsimulator, skip to next section.

.. note::
    The source instillation process differs for Macs with Intel processors and Macs with Apple
    silicon processors. Before continuing, check which processor type your computer has using
    `Apple's instructions <https://support.apple.com/en-us/HT211814>`__.

.. tabs::

  .. tab:: Linux
    :tabid: linus_source

    .. include:: source_install/linux.rst

  .. tab:: Mac OSX Intel
    :tabid: macosx_intel_source

    .. include:: source_install/macosx_intel.rst

  .. tab:: Mac OSX M1
    :tabid: macosx_m1_source

    .. include:: source_install/macosx_m1.rst

  .. tab:: Windows
    :tabid: windows_source

    .. include:: source_install/windows.rst

Build and install
"""""""""""""""""

Before building the development version of mrsimulator, install the development requirement
packages with pip. In the directory where your copy of mrsimulator is, type:

.. code-block:: bash

    $ pip install -r requirements-dev.txt
    $ pip install -e .

As always, if you get an error that you don’t have the permission to install the
package into the default site-packages directory, you may try installing by adding the
``--user`` option.


Note for the developers and contributors
''''''''''''''''''''''''''''''''''''''''

**Before commits**: Mrsimulator follows python community standards for writing code and
documentation. To help guide the developers and contributors towards these standards,
we have created a *.pre-commit-config.yaml* file, that when used with ``pre-commit``, will
inspect the code and document for issues.
Type ``pre-commit run`` before git commits to inspect the changes.

You can also set up the git hook script to automatically run *pre-commit* on git
commits with the ``pre-commit install``. Read more about
`pre-commit <https://pre-commit.com/#3-install-the-git-hook-scripts>`_.


**Running tests**: For unit tests, we use the pytest module. At the root directory
of the mrsimulator package folder, type

.. code-block:: bash

    $ pytest

which will run a series of tests.

**Building docs**: We use the sphinx python documentation generator for building docs.
Navigate to the *docs* folder within the mrsimulator package folder, and type,

.. code-block:: bash

    $ make html

The above command will build the documentation and store the build at
*mrsimulator/docs/_build/html*. Double click the *index.html* file within this
folder to view the offline documentation.

.. **Submitting pull requests** Make sure all the test pass and the documentation build
.. is successful before creating a pull request.

.. We recommend the
.. following C-compiler for the OS types:
.. - Mac OS - ``clang``
.. - Linux - ``gcc``
.. - Windows - ``msvc`` (https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2019)
