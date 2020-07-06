.. _czjzek_distribution:

Czjzek distribution
===================

Introduction
------------

A Czjzek distribution model [#f1]_ [#f2]_ is a random distribution of the second-rank
traceless symmetric tensors about a zero tensor. An explicit
form of a traceless symmetric second-rank tensor, :math:`{\bf S}`, in Cartesian basis,
follows,

.. math::
    {\bf S} = \left[
    \begin{array}{l l l}
    S_{xx} & S_{xy} & S_{xz} \\
    S_{xy} & S_{yy} & S_{yz} \\
    S_{xz} & S_{yz} & S_{zz}
    \end{array}
    \right],

where :math:`S_{xx} + S_{yy} + S_{zz} = 0`.
The elements of the above Cartesian tensor, :math:`S_{ij}`, can be decomposed into
second-rank irreducible spherical tensor components [#f3]_, :math:`R_{2,k}`, following

.. math::
    \begin{align}
    S_{xx} &= \frac{1}{2} (R_{2,2} + R_{2,-2}) - \frac{1}{\sqrt{6}} R_{2,0}, \\
    S_{xy} &= S_{yx} = -\frac{i}{2} (R_{2,2} - R_{2,-2}), \\
    S_{yy} &= -\frac{1}{2} (R_{2,2} + R_{2,-2}) - \frac{1}{\sqrt{6}} R_{2,0}, \\
    S_{xz} &= S_{zx} = -\frac{1}{2} (R_{2,1} - R_{2,-1}), \\
    S_{zz} &= \sqrt{\frac{2}{3}} R_{20}, \\
    S_{yz} &= S_{zy} = \frac{i}{2} (R_{2,1} + R_{2,-1}).
    \end{align}


In the Czjzek model, the distribution of the second-rank traceless symmetric tensor is
based on the assumption of a random distribution of the five irreducible spherical
tensor components, :math:`R_{2,k}`, drawn from an uncorrelated five-dimensional
multivariate normal distribution.
Since :math:`R_{2,k}` components are complex, random sampling is performed on the
equivalent real tensor components, which are a linear combination of :math:`R_{2,k}`,
and are given as

.. math::
    \begin{align}
    U_1 &= \frac{1}{\sqrt{6}} R_{2,0}, \\
    U_2 &= -\frac{1}{\sqrt{12}} (R_{2,1} - R_{2,-1}), \\
    U_3 &= \frac{i}{\sqrt{12}} (R_{2,1} + R_{2,-1}), \\
    U_4 &= -\frac{i}{\sqrt{12}} (R_{2,2} - R_{2,-2}), \\
    U_5 &= \frac{1}{\sqrt{12}} (R_{2,2} + R_{2,-2}),
    \end{align}

where :math:`U_i` forms an ortho-normal basis. The components, :math:`U_i`, are drawn
from a five-dimensional uncorrelated multivariate normal distribution with zero mean
and covariance matrix, :math:`\Lambda=\sigma^2 {\bf I}_5`, where :math:`{\bf I}_5` is a
:math:`5 \times 5` identity matrix and :math:`\sigma` is the standard deviation.

In terms of :math:`U_i`, the traceless second-rank symmetric Cartesian tensor elements,
:math:`S_{ij}`, follows

.. math::
    \begin{align}
    S_{xx} &= \sqrt{3} U_5 - U_1, \\
    S_{xy} &= S_{yx} = \sqrt{3} U_4, \\
    S_{yy} &= -\sqrt{3} U_5 - U_1, \\
    S_{xz} &= S_{zx} = \sqrt{3} U_2, \\
    S_{zz} &= 2 U_1, \\
    S_{yz} &= S_{zy} = \sqrt{3} U_3,
    \end{align}

and the explicit matrix form the :math:`{\bf S}` is

.. math::
    {\bf S} = \left[
    \begin{array}{l l l}
    \sqrt{3} U_5 - U_1   & \sqrt{3} U_4          & \sqrt{3} U_2 \\
    \sqrt{3} U_4         & -\sqrt{3} U_5 - U_1   & \sqrt{3} U_3 \\
    \sqrt{3} U_2         & \sqrt{3} U_3          & 2 U_1
    \end{array}
    \right].


Generating Czjzek distribution
------------------------------

Czjzek distribution of symmetric shielding tensors
''''''''''''''''''''''''''''''''''''''''''''''''''

To generate a Czjzek distribution, use the :func:`mrsimulator.models.czjzek_distribution`
function as follows,

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> from mrsimulator.models import czjzek_distribution
    >>> zeta_dist, eta_dist = czjzek_distribution(sigma=3.1415, n=10000)

The plot of the :math:`\zeta`-:math:`\eta` distribution.

.. plot::
    :format: python
    :context: close-figs
    :include-source:

    >>> import matplotlib.pyplot as plt # doctest: +SKIP
    >>> plt.scatter(zeta_dist, eta_dist, s=2) # doctest: +SKIP
    >>> plt.xlabel('$\zeta$ / ppm') # doctest: +SKIP
    >>> plt.ylabel('$\eta$') # doctest: +SKIP
    >>> plt.tight_layout() # doctest: +SKIP
    >>> plt.show() # doctest: +SKIP

.. The following is a simulation resulting from a Czjzek distribution of the symmetric
.. shielding tensors.

.. .. plot::
..     :format: doctest
..     :context: close-figs
..     :include-source:

..     >>> from mrsimulator import Simulator, SpinSystem, SpinSystem
..     >>> from mrsimulator.methods import BlochDecaySpectrum
..     ...
..     >>> systems = [
..     ...     SpinSystem(
..     ...         sites=[Site(isotope='13C', shielding_symmetric={'zeta': z, 'eta': e})],
..     ...         abundance=1e-4,
..     ...     )
..     ...     for z, e in zip(zeta_dist, eta_dist)
..     ... ]
..     >>> method = BlochDecaySpectrum(
..     ...     channels=['13C'],
..     ...     spectral_dimensions=[{'count': 1024, 'spectral_width': 1e4}]
..     ... )
..     ...
..     >>> sim = Simulator()
..     >>> sim.spin_systems = systems
..     >>> sim.methods = [method]
..     >>> sim.run()
..     ...
..     >>> ax = plt.gca(projection='csdm') # doctest: +SKIP
..     >>> ax.plot(sim.methods[0].simulation) # doctest: +SKIP
..     >>> ax.invert_xaxis() # doctest: +SKIP
..     >>> plt.tight_layout() # doctest: +SKIP
..     >>> plt.show() # doctest: +SKIP


----


Czjzek distribution of symmetric quadrupolar tensors
''''''''''''''''''''''''''''''''''''''''''''''''''''

The Czjzek distribution of symmetric quadrupolar tensors follows a similar setup as the
Czjzek distribution of symmetric shielding tensors.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> from mrsimulator.models import czjzek_distribution
    >>> Cq_dist, eta_dist = czjzek_distribution(sigma=1.12, n=10000)

The plot of the :math:`\zeta`-:math:`\eta` distribution.

.. plot::
    :format: python
    :context: close-figs
    :include-source:

    >>> import matplotlib.pyplot as plt # doctest: +SKIP
    >>> plt.scatter(Cq_dist, eta_dist, s=2) # doctest: +SKIP
    >>> plt.xlabel('$C_q$ / MHz') # doctest: +SKIP
    >>> plt.ylabel('$\eta$') # doctest: +SKIP
    >>> plt.tight_layout() # doctest: +SKIP
    >>> plt.show() # doctest: +SKIP

.. The following is a simulation resulting from a Czjzek distribution of the quadrupolar
.. tensors.

.. .. plot::
..     :format: doctest
..     :context: close-figs
..     :include-source:

..     >>> from mrsimulator import Simulator, SpinSystem, SpinSystem
..     >>> from mrsimulator.methods import BlochDecayCentralTransitionSpectrum
..     ...
..     >>> systems = [
..     ...     SpinSystem(
..     ...         sites=[Site(isotope='71Ga', quadrupolar={'Cq': c * 1e6, 'eta': e})],
..     ...         abundance=1e-4,
..     ...     )
..     ...     for c, e in zip(Cq_dist, eta_dist)
..     ... ]
..     >>> method = BlochDecayCentralTransitionSpectrum(
..     ...     channels=['71Ga'],
..     ...     magnnetic_flux_density=4.6,
..     ...     spectral_dimensions=[
..     ...         {'count': 1024, 'spectral_width': 1.2e5, 'reference_offset': -1e4}
..     ...     ]
..     ... )
..     ...
..     >>> sim = Simulator()
..     >>> sim.spin_systems = systems
..     >>> sim.methods = [method]
..     >>> sim.run()
..     ...
..     >>> ax = plt.gca(projection='csdm') # doctest: +SKIP
..     >>> ax.plot(sim.methods[0].simulation) # doctest: +SKIP
..     >>> ax.invert_xaxis() # doctest: +SKIP
..     >>> plt.tight_layout() # doctest: +SKIP
..     >>> plt.show() # doctest: +SKIP

References
----------

.. [#f1] Czjzek, G., Fink, J., Götz, F., Schmidt, H., Coey, J. M. D., Atomic
    coordination and the distribution of electric field gradients in amorphous solids
    Phys. Rev. B (1981) **23** 2513-30.
    `DOI: 10.1103/PhysRevB.23.2513 <https://doi.org/10.1103/PhysRevB.23.2513>`_

.. [#f2] Caër, G.L., Bureau, B., Massiot, D., An extension of the Czjzek model for the
    distributions of electric field gradients in disordered solids and an application
    to NMR spectra of 71Ga in chalcogenide glasses. Journal of Physics: Condensed
    Matter, (2010), **22**.
    `DOI: 10.1088/0953-8984/22/6/065402 <https://doi.org/10.1088/0953-8984/22/6/065402>`_

.. [#f3] Grandinetti, P. J., Ash, J. T., Trease, N. M. Symmetry pathways in solid-state
    NMR, PNMRS 2011 **59**, *2*, 121-196.
    `DOI: 10.1016/j.pnmrs.2010.11.003 <https://doi.org/10.1016/j.pnmrs.2010.11.003>`_


.. minigallery:: mrsimulator.models.czjzek_distribution
    :add-heading: Mini-gallery using czjzek distributions
    :heading-level: -
