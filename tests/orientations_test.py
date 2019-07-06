import numpy as np
from mrsimulator.python.orientation import trig_of_polar_angles_and_amplitudes
from mrsimulator.python.orientation import triangle_interpolation

import mrsimulator.sandbox as clib


def test_polar_coordinates_and_amps():
    nt = 64
    cos_alpha_py, cos_beta_py, amp_py = trig_of_polar_angles_and_amplitudes(nt)
    cos_alpha_c, cos_beta_c, amp_c = clib.trig_of_polar_angles_and_amplitudes(nt)

    assert np.allclose(cos_alpha_py, cos_alpha_c, atol=1e-15)
    assert np.allclose(cos_beta_py, cos_beta_c, atol=1e-15)
    assert np.allclose(amp_py, amp_c, atol=1e-15)


def test_triangle_interpolation():
    f = np.asarray([10, 30, 50])
    amp_py = np.zeros(100)
    triangle_interpolation(f, amp_py)

    amp_c = np.zeros(100)
    clib.triangle_interpolation(f, amp_c)

    assert np.allclose(amp_py, amp_c, atol=1e-15)
