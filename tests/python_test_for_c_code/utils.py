# -*- coding: utf-8 -*-
import numpy as np

__author__ = "Deepansh J. Srivastava"
__email__ = "deepansh2012@gmail.com"


def pre_phase_components(number_of_sidebands, spin_frequency):
    r"""
    Calculates a 2D matrix `pre_phase` of shape (number_of_sidebands, 9) with
    the following definition

    `pre_phase = $2i\pi \frac{(\exp(i m \omega_r t) - 1)}{i m \omega_r}$`

    where m goes from -4 to 4 and t is the time points of size
    `number_of_points`. The time increment is given as
    `tau = 1.0 / (number_of_sidebands * spin_frequency)`.
    """
    n = int(number_of_sidebands)
    m_wr = spin_frequency * 2.0 * np.pi * (np.arange(9) - 4.0)
    time = np.arange(n) / (number_of_sidebands * spin_frequency)
    time, m_wr = np.meshgrid(time, m_wr)
    pht = 1j * time * m_wr
    pre_phase = np.empty(pht.shape, np.complex128)
    np.expm1(pht, out=pre_phase)
    pre_phase[0:4] *= (2.0 * np.pi) / m_wr[0:4]
    pre_phase[5:9] *= (2.0 * np.pi) / m_wr[5:9]
    pre_phase[4] = 0.0
    return pre_phase
