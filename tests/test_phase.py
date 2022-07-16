"""Test for c functions."""
import mrsimulator.tests.tests as clib
import numpy as np

from .python_test_for_c_code.utils import pre_phase_components


def test_phase_components_1():
    number_of_sidebands = 64
    spin_frequency = 10
    pre_phase_py = pre_phase_components(number_of_sidebands, spin_frequency)
    pre_phase_c = clib.pre_phase_components(number_of_sidebands, spin_frequency)

    assert np.allclose(pre_phase_c, 2 * pre_phase_py[:4, :])


def test_phase_components_2():
    number_of_sidebands = 32
    spin_frequency = 10e3
    pre_phase_py = pre_phase_components(number_of_sidebands, spin_frequency)
    pre_phase_c = clib.pre_phase_components(number_of_sidebands, spin_frequency)

    assert np.allclose(pre_phase_c, 2 * pre_phase_py[:4, :])


def test_phase_components_3():
    number_of_sidebands = 128
    spin_frequency = 12.5e3
    pre_phase_py = pre_phase_components(number_of_sidebands, spin_frequency)
    pre_phase_c = clib.pre_phase_components(number_of_sidebands, spin_frequency)

    assert np.allclose(pre_phase_c, 2 * pre_phase_py[:4, :])
