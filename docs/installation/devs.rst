For contributors and developers
===============================

Make your own copy of mrsimulator on GitHub
'''''''''''''''''''''''''''''''''''''''''''

Making a copy of someone's code on GitHub is the same as making a *fork*. A fork is a
complete copy of the code and all of its revision history.

1. Log into a `GitHub account <https://github.com>`_.
2. Go to the `mrsimulator Github <https://github.com/DeepanshS/mrsimulator>`_ home page.
3. Click on the *fork* button.

You will see a short animation of Octocat scanning a book on a flatbed scanner. After
that, you should find yourself at the home page for your forked copy of mrsimulator.


Create a development environment
''''''''''''''''''''''''''''''''

It is good practice to create separate virtual python environments for packages when
in developer mode. The following is an example of a Conda environment.

.. code-block:: bash

    $ conda create -n mrsimulator-dev python=3.7

The above command will create a new python3.7 environment named *mrsimulator-dev*. To
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
    git remote add mrsimulator git://github.com/DeepanshS/mrsimulator.git

You can check that everything is set up properly so far by asking git to show you all of the
remotes it knows about for your local repository of mrsimulator with ``git remote -v``, which
should display something like

.. code-block:: bash

    mrsimulator   git://github.com/DeepanshS/mrsimulator.git (fetch)
    mrsimulator   git://github.com/DeepanshS/mrsimulator.git (push)
    origin     git@github.com:your-user-name/mrsimulator.git (fetch)
    origin     git@github.com:your-user-name/mrsimulator.git (push)


Build the development version of mrsimulator
''''''''''''''''''''''''''''''''''''''''''''

OS-dependent prerequisites
""""""""""""""""""""""""""

.. note::
    Installing OS-dependent prerequisites is a one-time process. If you are
    upgrading to a newer version of mrsimulator, skip to next section.

.. tabs::

  .. tab:: Linux
    :tabid: linus_source

    .. include:: source_install/linux.rst

  .. tab:: Mac OSX
    :tabid: macosx_source

    .. include:: source_install/macosx.rst

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
