..  _requirements:

Package dependencies
====================

``mrsimulator`` is compatible with the following operating systems:

- MacOS 10.15 or later
- Windows 7 or later
- Most releases of Linux

``mrsimulator`` depends on the following packages:

**Required packages**

- `NumPy>=1.17 <http://www.numpy.org>`_
- openblas
- cython>=0.29.14
- typing-extensions>=3.7
- `matplotlib>=3.4 <https://matplotlib.org>`_ for figures and visualization,
- monty>=2.0.4
- `csdmpy>=0.4.1 <https://csdmpy.readthedocs.io/en/stable/>`_
- `pydantic>=1.9 <https://pydantic-docs.helpmanual.io>`_

**Other packages**

- pytest>=4.5.0 for unit tests.
- pre-commit for code formatting
- sphinx>=2.0 for generating the documentation
- sphinxjp.themes.basicstrap for documentation.
- breathe==4.26 for generating C documentation
- sphinx-copybutton
