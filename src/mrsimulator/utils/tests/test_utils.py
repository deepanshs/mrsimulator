# -*- coding: utf-8 -*-
import csdmpy as cp
import numpy as np
from mrsimulator.utils import get_spectral_dimensions


def test_get_spectral_dimensions():
    # 1
    csdm = cp.as_csdm(np.arange(10))
    csdm.x[0] = cp.LinearDimension(count=10, increment="1 Hz", complex_fft=True)

    res = {"count": 10, "spectral_width": 10, "reference_offset": 0}
    assert get_spectral_dimensions(csdm)[0] == res

    # 2
    csdm.x[0] = cp.LinearDimension(
        count=10, increment="1 Hz", coordinates_offset="-0.05 kHz", complex_fft=True
    )

    res = {"count": 10, "spectral_width": 10, "reference_offset": -50}
    assert get_spectral_dimensions(csdm)[0] == res

    # 3
    csdm.x[0] = cp.LinearDimension(count=10, increment="1 Hz")

    res = {"count": 10, "spectral_width": 10, "reference_offset": 5}
    assert get_spectral_dimensions(csdm)[0] == res

    # 4
    csdm.x[0] = cp.LinearDimension(count=10, increment="1 Hz", label="1H")

    res = {"count": 10, "spectral_width": 10, "reference_offset": 5, "label": "1H"}
    assert get_spectral_dimensions(csdm)[0] == res

    # 5
    csdm.x[0] = cp.LinearDimension(count=10, increment="-1 Hz", label="1H")

    res = {"count": 10, "spectral_width": 10, "reference_offset": -4, "label": "1H"}
    assert get_spectral_dimensions(csdm)[0] == res

    # 6
    csdm.x[0] = cp.LinearDimension(
        count=10, increment="-1 Hz", coordinates_offset="-10 Hz", label="1H"
    )

    res = {"count": 10, "spectral_width": 10, "reference_offset": -14, "label": "1H"}
    assert get_spectral_dimensions(csdm)[0] == res

    # 7
    csdm.x[0] = cp.LinearDimension(
        count=10,
        increment="-1 Hz",
        coordinates_offset="-10 Hz",
        origin_offset="100 MHz",
        label="1H",
    )

    res = {
        "count": 10,
        "spectral_width": 10,
        "reference_offset": -14,
        "origin_offset": 100e6,
        "label": "1H",
    }
    assert get_spectral_dimensions(csdm)[0] == res
