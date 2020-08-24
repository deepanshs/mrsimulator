.. _czjzek_distribution:

Czjzek distribution
-------------------

A Czjzek distribution model is a random distribution of the second-rank traceless
symmetric tensors about a zero tensor. See :ref:`czjzek_model` and references within
for a brief description of the model.

Czjzek distribution of symmetric shielding tensors
''''''''''''''''''''''''''''''''''''''''''''''''''

To generate a Czjzek distribution, use the :class:`~mrsimulator.models.CzjzekDistribution`
class as follows.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> from mrsimulator.models import CzjzekDistribution
    >>> cz_model = CzjzekDistribution(sigma=0.8)

The `CzjzekDistribution` class accepts a single argument, `sigma`, which is the standard
deviation of the second-rank traceless symmetric tensor parameters. In the above example,
we create ``cz_model`` as an instance of the CzjzekDistribution class with
:math:`\sigma=0.8`.

Note, `cz_model` is only a class instance of the Czjzek distribution. You can either
draw random points from this distribution or generate a probability distribution
function. Let's first draw points from this distribution, using the
:meth:`~mrsimulator.models.CzjzekDistribution.rvs` method of the instance.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> zeta_dist, eta_dist = cz_model.rvs(n=50000)

In the above example, we draw `n=50000` random points of the distribution. The output
``zeta_dist`` and ``eta_dist`` hold the coordinates of the points.
The scatter plot of these coordinates is shown below.

.. plot::
    :format: python
    :context: close-figs
    :include-source:

    >>> import matplotlib.pyplot as plt # doctest: +SKIP
    >>> plt.scatter(zeta_dist, eta_dist, s=4, alpha=0.02) # doctest: +SKIP
    >>> plt.xlabel('$\zeta$ / ppm') # doctest: +SKIP
    >>> plt.ylabel('$\eta$') # doctest: +SKIP
    >>> plt.xlim(-15, 15) # doctest: +SKIP
    >>> plt.ylim(0, 1) # doctest: +SKIP
    >>> plt.tight_layout() # doctest: +SKIP
    >>> plt.show() # doctest: +SKIP

----

Czjzek distribution of symmetric quadrupolar tensors
''''''''''''''''''''''''''''''''''''''''''''''''''''

The Czjzek distribution of symmetric quadrupolar tensors follows a similar setup as the
Czjzek distribution of symmetric shielding tensors, except we assign the outputs to Cq
and :math:`\eta_q`. In the following example, we generate the probability distribution
function using the :meth:`~mrsimulator.models.CzjzekDistribution.pdf` method.

.. plot::
    :format: doctest
    :context: close-figs
    :include-source:

    >>> import numpy as np
    >>> Cq_range = np.arange(100)*0.3 - 15 # pre-defined Cq range in MHz.
    >>> eta_range = np.arange(21)/20  # pre-defined eta range.
    >>> Cq, eta, amp = cz_model.pdf(pos=[Cq_range, eta_range])

To generate a probability distribution, we need to define a grid system over which the
distribution probabilities will be evaluated. We do so by defining the range of coordinates
along the two dimensions. In the above example, ``Cq_range`` and ``eta_range`` are the
range of :math:`\text{Cq}` and :math:`\eta_q` coordinates, which is then given as the
argument to the :meth:`~mrsimulator.models.CzjzekDistribution.pdf` method. The output
``Cq``, ``eta``, and ``amp`` hold the two coordinates and amplitude, respectively.

The plot of the Czjzek probability distribution is shown below.

.. plot::
    :format: python
    :context: close-figs
    :include-source:

    >>> import matplotlib.pyplot as plt # doctest: +SKIP
    >>> plt.contourf(Cq, eta, amp, levels=10) # doctest: +SKIP
    >>> plt.xlabel('$C_q$ / MHz') # doctest: +SKIP
    >>> plt.ylabel('$\eta$') # doctest: +SKIP
    >>> plt.tight_layout() # doctest: +SKIP
    >>> plt.show() # doctest: +SKIP

.. note::
    The ``pdf`` method of the instance generates the probability distribution function
    by first drawing random points from the distribution and then binning it
    onto a pre-defined grid.

.. minigallery:: mrsimulator.models.CzjzekDistribution
    :add-heading: Mini-gallery using the Czjzek distributions
    :heading-level: '
