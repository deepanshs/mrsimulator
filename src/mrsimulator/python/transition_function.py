# -*- coding: utf-8 -*-
"""Transition symmetry functions from irreducible tensors.

.. seealso::

    `Grandinetti et. al. <https://doi.org/10.1016/j.pnmrs.2010.11.003>`_
    **Symmetry pathways in solid-state NMR.**
"""

__all__ = ["p", "d", "dIS", "f", "quad_ci"]


def p(mf, mi):
    r"""
    Transition symmetry function from irreducible 1st rank tensor.
    The function evaluates,

    .. math::

        \mathbb{p}(m_f, m_i) &= \left< m_f | T_{10} | m_f \right> -
                                \left< m_i | T_{10} | m_i \right> \\
        &= m_f - m_i

    where :math:`T_{10}` is the irreducible 1st rank tensor operator in the
    rotating tilted frame.

    :param float mf: The quantum number associated with the final energy state.
    :param float mi: The quantum number associated with the initial energy
                     state.
    :rtype: float
    :return: The transition symmetry function p.
    """
    return mf - mi


def d(mf, mi):
    r"""
    Transition symmetry function from irreducible 2nd rank tensor.
    The function evaluates,

    .. math::

        \mathbb{d}(m_f, m_i) &= \left< m_f | T_{20} | m_f \right> -
                                \left< m_i | T_{20} | m_i \right> \\
        &= \sqrt{\frac{3}{2}} \left(m_f^2 - m_i^2 \right)

    where :math:`T_{20}` is the irreducible 2nd rank tensor operator in the
    rotating tilted frame.

    :param float mf: The quantum number associated with the final energy state.
    :param float mi: The quantum number associated with the initial energy
                     state.
    :rtype: float
    :return: The transition symmetry function d.
    """
    return 1.2247448714 * (mf * mf - mi * mi)


def f(mf, mi, spin):
    r"""
    Transition symmetry function from irreducible 3rd rank tensor.
    The function evaluates,

    .. math::

        \mathbb{f}(m_f, m_i) &= \left< m_f | T_{30} | m_f \right> -
                                \left< m_i | T_{30} | m_i \right> \\
        &= \frac{1}{\sqrt{10}} [5(m_f^3 - m_i^3) + (1 - 3I(I+1))(m_f-m_i)]

    where :math:`T_{30}` is the irreducible 3rd rank tensor operator in the
    rotating tilted frame.

    :param float mf: The quantum number associated with the final energy state.
    :param float mi: The quantum number associated with the initial energy
                     state.
    :rtype: float
    :return: The transition symmetry function f.
    """
    f_ = 1.0 - 3.0 * spin * (spin + 1.0)
    f_ *= mf - mi
    f_ += 5.0 * (mf * mf * mf - mi * mi * mi)
    f_ *= 0.316227766
    return f_


def quad_ci(mf, mi, spin):
    f_ = f(mf, mi, spin)
    p_ = p(mf, mi)
    temp = spin * (spin + 1.0) - 0.75
    c0 = 0.3577708764 * temp * p_ + 0.8485281374 * f_
    c2 = 0.1069044968 * temp * p_ + -1.0141851057 * f_
    c4 = -0.1434274331 * temp * p_ + -1.2850792082 * f_
    return c0, c2, c4


def dIS(mIf, mIi, mSf, mSi):
    return mIf * mSf - mIi * mSi
