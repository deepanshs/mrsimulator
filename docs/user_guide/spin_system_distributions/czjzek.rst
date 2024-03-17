.. _czjzek_distribution:

Czjzek distribution
-------------------

The Czjzek distribution models random variations of second-rank traceless
symmetric tensors about zero, i.e., a tensor with zeta of zero. An analytical expression
for the Czjzek distribution exists (cite) which follows

.. math::
    f(\zeta, \eta, \sigma) = \eta \left(1-\frac{\eta^2}{9}\right)\frac{\zeta^4}{32\sigma^5 \sqrt{2 \pi}} \times \exp\left(-\frac{\zeta^2}{8\sigma^2}\left(1+\frac{\eta^2}{3}\right)\right),

where :math:`\zeta` and :math:`\eta` are the Haberlen components of the tensor and :math:`\sigma` is the Czjzek width parameter. See :ref:`czjzek_model` for a further mathematical description of the model.

The remainder of this page quickly describes how to generate Czjzek distributions and generate
:py:class:`~mrsimulator.spin_system.SpinSystem` objects from these distributions. Also, look at the
gallery examples using the Czjzek distribution listed at the bottom of this page.

Creating and sampling a Czjzek distribution
'''''''''''''''''''''''''''''''''''''''''''

To generate a Czjzek distribution, use the :py:class:`~mrsimulator.models.CzjzekDistribution`
class as follows.

.. plot::
    :context: reset

    from mrsimulator.models import CzjzekDistribution

    cz_model = CzjzekDistribution(sigma=0.8)

The **CzjzekDistribution** class accepts the argument, ``sigma``, which is the standard
deviation of the second-rank traceless symmetric tensor parameters. In the above example,
we create ``cz_model`` as an instance of the CzjzekDistribution class with
:math:`\sigma=0.8`.

Note, ``cz_model`` is only a class instance of the Czjzek distribution. You can either
draw random points from this distribution or generate a probability distribution
function. Let's first draw points from this distribution, using the
:py:meth:`~mrsimulator.models.CzjzekDistribution.rvs` method of the instance.

.. plot::
    :context: close-figs

    zeta_dist, eta_dist = cz_model.rvs(size=50000)

In the above example, we draw *50000* random points of the distribution. The output ``zeta_dist`` and ``eta_dist`` hold the tensor parameter coordinates of the points, defined in the Haeberlen convention. It is further assumed that the points in ``zeta_dist`` are in units of ``ppm`` while ``eta_dist`` has values since :math:`\eta` is dimensionless. The scatter plot of these coordinates is shown below.

.. skip: next

.. plot::
    :context: close-figs
    :caption: Random sampling of Czjzek distribution of shielding tensors.

    import matplotlib.pyplot as plt

    plt.scatter(zeta_dist, eta_dist, s=4, alpha=0.02)
    plt.xlabel("$\zeta$ / ppm")
    plt.ylabel("$\eta$")
    plt.xlim(-15, 15)
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.show()

Creating and sampling a Czjzek distribution in polar coordinates
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

The :py:class:`~mrsimulator.models.czjzek.CzjzekDistribution` class also supports sampling tensors in polar coordinates. The logic behind transforming from a :math:`\zeta`-:math:`\eta` Cartesian grid is further described in mrinversion (cite), and the following definitions are used

.. math::

    \begin{split}r_\zeta = \left| \zeta \right| ~~~~\text{and}~~~~
    \theta = \left\{ \begin{array}{l r}
                \frac{\pi}{4} \eta      &: \zeta \le 0, \\
                \frac{\pi}{2} \left(1 - \frac{\eta}{2} \right) &: \zeta > 0.
            \end{array}
            \right.\end{split}

Because Cartesian grids are more manageable in computation, the above polar piece-wise grid is re-express as the x-y Cartesian grid following,

.. math::

    x = r_\zeta \cos\theta ~~~~\text{and}~~~~ y = r_\zeta \sin\theta.

Below, we create another instance of the :py:class:`~mrsimulator.models.czjzek.CzjzekDistribution`
class with the same value of :math:`sigma=0.8`, but we now also include the argument ``polar=True``
which means the :py:meth:`~mrsimulator.models.CzjzekDistribution.rvs` will sample x and y points.

.. skip: next

.. plot::
    :context: close-figs
    :caption: Random sampling of Czjzek distribution of shielding tensors in polar coordinates.

    cz_model_polar = CzjzekDistribution(sigma=0.8, polar=True)

    # Sample (x, y) points
    x_dist, y_dist = cz_model_polar.rvs(size=50000)

    # Plot the distribution
    plt.figure(figsize=(4, 4))
    plt.scatter(x_dist, y_dist, s=4, alpha=0.02)
    plt.xlabel("$x$ / ppm")
    plt.ylabel("$y$ / ppm")
    plt.xlim(0, 8)
    plt.ylim(0, 8)
    plt.tight_layout()
    plt.show()


----

Generating probability distribution functions from a Czjzek model
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

The :py:meth:`~mrsimulator.models.CzjzekDistribution.pdf` instance method will generate a
probability distribution function on the supplied grid using the analytical function defined above.
The provided grid -- passed to the ``pos`` keyword argument -- needs to be defined in either
Cartesian or polar coordinates depending on if the
:py:attr:`~mrsimulator.models.CzjzekDistribution.polar` attribute is ``True`` or ``False``.

Below, we generate and plot a probability distribution on a :math:`\zeta`-:math:`\eta` Cartesian
grid where ``zeta_range`` and ``eta_range`` define the desired coordinates in each dimension of the
grid system.

.. plot::
    :context: close-figs

    import numpy as np

    cz_model = CzjzekDistribution(sigma=1.2, polar=False)  # sample in (zeta, eta)

    zeta_range = np.linspace(-12, 12, num=200)  # pre-defined zeta range in units of ppm
    eta_range = np.linspace(0, 1, num=50)       # pre-defined eta range.
    zeta_grid, eta_grid, amp = cz_model.pdf(pos=[zeta_range, eta_range])

Here, ``zeta_grid`` and ``eta_grid`` are numpy arrays defining a set of pair-wise points on the
grid system, and ``amp`` is another numpy array holding the probability density at each point
on the grid. Below, the distribution is plotted

.. skip: next

.. plot::
    :context: close-figs
    :caption: Czjzek Distribution of shielding tensors.


    plt.contourf(zeta_grid, eta_grid, amp, levels=10)
    plt.xlabel("$\zeta$ / ppm")
    plt.ylabel("$\eta$")
    plt.tight_layout()
    plt.show()

---

The probability distribution function can also be generated in polar coordinates. The workflow
is the same, except we now define an (x, y) grid system using the variables ``x_range``
and ``y_range``. The code to generate and plot the polar Czjzek distribution is shown below.

.. skip: next

.. plot::
    :context: close-figs
    :caption: Czjzek Distribution of shielding tensors in polar coordinates.

    cz_model_polar = CzjzekDistribution(sigma=1.2, polar=True)  # sample in (x, y)

    x_range = np.linspace(0, 10, num=150)
    y_range = np.linspace(0, 10, num=150)
    x_grid, y_grid, amp = cz_model_polar.pdf(pos=[x_range, y_range])

    plt.figure(figsize=(4, 4))
    plt.contourf(x_grid, y_grid, amp, levels=10)
    plt.xlabel("$x$ / ppm")
    plt.ylabel("$y$ / ppm")
    plt.tight_layout()
    plt.show()


Distributions of shielding and quadrupolar tensors and a note on units
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

The :py:class:`~mrsimulator.models.CzjzekDistribution` class can be used to generate
distributions for both symmetric chemical shielding tensors and electric field gradient
tensors. It is important to note the Czjzek model is only aware of the Haberlen components
of the tensors and not the units of the tensor. In the above code cells, we generated
distributions for symmetric shielding tensors and assumed all units for :math:`\zeta` were
in ppm.

Quadrupolar tensors, defined using values of :math:`C_q` in MHz and unitless :math:`\eta`,
can also be drawn from the Czjzek distribution in the same manner; however, the dimensions
are assumed to be in units of MHz. The following code draws a distribution of quadrupolar
tensor parameters.

.. skip: next

.. plot::
    :context: close-figs

    Cq_range = np.linspace(-12, 12, num=200)  # pre-defined Cq range in units of MHz
    eta_range = np.linspace(0, 1, num=50)     # pre-defined eta range.
    Cq_grid, eta_grid, amp = cz_model.pdf(pos=[Cq_range, eta_range])


the units for ``Cq_range`` and ``Cq_grid`` are assumed in MHz. Similarly, x and y are assumed to
be in units of MHz when sampling quadrupolar tensors in polar coordinates.

.. skip: next

.. plot::
    :context: close-figs

    x_range = np.linspace(0, 10, num=150)  # pre-defined x grid in units of MHz
    y_range = np.linspace(0, 10, num=150)  # pre-defined y grid in units of MHz
    x_grid, y_grid, amp = cz_model_polar.pdf(pos=[x_range, y_range])

Generating a list of SpinSystem objects from a Czjzek model
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

The utility function :py:meth:`~mrsimulator.utils.collection.single_site_system_generator`, further
described in :ref:`single_site_system_generator_documentation`, can be used in conjunction with
the :py:class:`~mrsimulator.models.CzjzekDistribution` class to generate a set of spin systems whose
tensor parameters follow the Czjzek distribution.

.. plot::
    :context: close-figs

    from mrsimulator.utils.collection import single_site_system_generator

    # Distribution of quadrupolar tensors
    cz_model = CzjzekDistribution(sigma=0.7)
    Cq_range = np.linspace(0, 10, num=100)
    eta_range = np.linspace(0, 1, num=50)

    # Create (Cq, eta) grid points and amplitude
    Cq_grid, eta_grid, amp = cz_model.pdf(pos=[Cq_range, eta_range])

    sys = single_site_system_generator(
        isotope="27Al",
        quadrupolar={"Cq": Cq_grid * 1e6, "eta": eta_grid},  # Cq argument in units of Hz
        abundance=amp,
    )

A spin system will be generated for each point on the :math:`\zeta`-:math:`\eta` grid, and the
abundance of each spin system matches the amplitude of the Czjzek distribution. When working in
polar coordinates, the set of :math:`\left(x, y\right)` coordinates needs to be transformed into
a set of :math:`\left(\zeta, \eta\right)` coordinates before being passed to the
:py:meth:`~mrsimulator.utils.collection.single_site_system_generator` function. The utility
function :py:meth:`~mrsimulator.utils.x_y_to_zeta_eta` performs this transformation, as shown
below.

.. plot::
    :context: close-figs

    from mrsimulator.models.utils import x_y_to_zeta_eta

    # Sample distribution of shielding tensors in polar coords
    cz_model_polar = CzjzekDistribution(sigma=0.7, polar=True)
    x_range = np.linspace(0, 10, num=150)
    y_range = np.linspace(0, 10, num=150)

    # Create (x, y) grid points and amplitude
    x_grid, y_grid, amp = cz_model_polar.pdf(pos=[x_range, y_range])

    # To transformation (x, y) -> (zeta, eta)
    zeta_grid, eta_grid = x_y_to_zeta_eta(x_grid, y_grid)

---

.. minigallery:: mrsimulator.models.CzjzekDistribution
    :add-heading: Mini-gallery using the Czjzek distributions
    :heading-level: '
