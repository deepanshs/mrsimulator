import mrsimulator.sandbox as clib
from mrsimulator.python.utils import pre_phase_components
import numpy as np


def test_phase_components():
    number_of_sidebands = 64
    spin_frequency = 10
    pre_phase_py = pre_phase_components(number_of_sidebands, spin_frequency)
    pre_phase_c = clib.pre_phase_components(number_of_sidebands, spin_frequency)

    assert np.allclose(pre_phase_c, pre_phase_py)
