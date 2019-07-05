import numpy as np


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
    pht = time * m_wr
    scale = (2.0 * np.pi) / m_wr
    pre_phase = scale * (np.cos(pht) + 1j * np.sin(pht) - 1.0)
    pre_phase[4] = 0.0
    return pre_phase
