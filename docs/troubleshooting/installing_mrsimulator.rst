.. _installing_mrsimulator_troubleshooting

Errors and Issues when installing ``mrsimulator``
-------------------------------------------------

Python 3.10 incompatibility
"""""""""""""""""""""""""""

Some libraries ``mrsimulator`` depends on are not yet compatible with ``Python 3.10``. Running

.. code-block:: shell

    $ pip install mrsimulator

on ``Python 3.10`` will get stuck in a loop and the following error will appear:

.. code-block:: shell

    ERROR: Command errored out with exit status 1:

    ---  more lines omitted ---

    numpy.distutils.system_info.NotFoundError: No BLAS/LAPACK libraries found. Note: Accelerate is no longer supported.
    To build Scipy from sources, BLAS & LAPACK libraries need to be installed.
    See site.cfg.example in the Scipy source directory and
    https://docs.scipy.org/doc/scipy/reference/building/index.html for details.

After stopping the installation process, check the version of python by running

.. code-block:: shell

    $ python -V

in the command line. If ``Python 3.10.x`` is printed, please see details for `setting up a
environment manager <_package_manager_troubleshooting>`__ to create an environment with a lower
Python version.

Have a different installation issue?
""""""""""""""""""""""""""""""""""""

If you encountered a different encountered a different issue when installing ``mrsimulator``,
please submit an issue our `Github issue tracker <https://github.com/DeepanshS/mrsimulator/issues>`_.

Github is the preferred way for reporting, but issues can also be reported by
directly contacting `Matthew Giammar <mailto:giammar.7@osu.edu>`__.
