
.. _models:

Models
======

.. _czjzek_model:

Czjzek distribution
-------------------

A Czjzek distribution model [#f1]_ is a random distribution of the second-rank
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

and the explicit matrix form of :math:`{\bf S}` is

.. math::
    {\bf S} = \left[
    \begin{array}{l l l}
    \sqrt{3} U_5 - U_1   & \sqrt{3} U_4          & \sqrt{3} U_2 \\
    \sqrt{3} U_4         & -\sqrt{3} U_5 - U_1   & \sqrt{3} U_3 \\
    \sqrt{3} U_2         & \sqrt{3} U_3          & 2 U_1
    \end{array}
    \right].

In a shorthand notation, we denote a Czjzek distribution of second-rank traceless
symmetric tensor as :math:`S_C(\sigma)`.

.. _ext_czjzek_model:

Extended Czjzek distribution
----------------------------

An Extended Czjzek distribution model [#f2]_ is a random perturbation of the second-rank
traceless symmetric tensors about a non-zero tensor, which is given as

.. math::
    S_T = S(0) + \rho S_C(\sigma=1),

where :math:`S_T` is the total tensor, :math:`S(0)` is the non-zero dominant second-rank
tensor, :math:`S_C(\sigma=1)` is the Czjzek random model attributing to the random
perturbation of the tensor about the dominant tensor, :math:`S(0)`, and :math:`\rho` is
the size of the perturbation. Note, in the above equation, the :math:`\sigma` parameter
from the Czjzek random model, :math:`S_C`, has no meaning and is set to one. The factor,
:math:`\rho`, is defined as

.. math::
    \rho = \frac{||S(0)|| \epsilon}{\sqrt{30}},

where :math:`\|S(0)\|` is the 2-norm of the dominant tensor, and :math:`\epsilon` is a
fraction.

----

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
