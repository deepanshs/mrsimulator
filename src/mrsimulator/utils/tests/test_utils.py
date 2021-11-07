# -*- coding: utf-8 -*-
import csdmpy as cp
import numpy as np
from mrsimulator.method import SpectralDimension
from mrsimulator.utils import get_spectral_dimensions


def assertion(csdm, res, test_ref=""):
    assert get_spectral_dimensions(csdm)[0] == res, f"failed ref {test_ref}"

    csdm_coordinates = csdm.x[0].coordinates.to("Hz").value
    if csdm.x[0].increment < 0:
        csdm_coordinates = csdm_coordinates[::-1]
    spectral_dimension = SpectralDimension(**res)
    mrsim_coordinates = spectral_dimension.coordinates_Hz()
    assert np.allclose(csdm_coordinates, mrsim_coordinates), f"failed ref {test_ref}"


def test_get_spectral_dimensions_even_count_positive_increment():
    # even # complex_fft
    csdm = cp.as_csdm(np.arange(10))
    csdm.x[0] = cp.LinearDimension(count=10, increment="1 Hz", complex_fft=True)
    res = {"count": 10, "spectral_width": 10, "reference_offset": 0}
    assertion(csdm, res, "even complex_fft")

    csdm.x[0] = cp.LinearDimension(count=10, increment="3.4 Hz", complex_fft=True)
    res = {"count": 10, "spectral_width": 34, "reference_offset": 0}
    assertion(csdm, res, "even complex_fft")

    # even # complex_fft # coordinates_offset
    csdm.x[0] = cp.LinearDimension(
        count=10, increment="1 Hz", coordinates_offset="-0.05 kHz", complex_fft=True
    )
    res = {"count": 10, "spectral_width": 10, "reference_offset": -50}
    assertion(csdm, res, "even complex_fft coordinates_offset")

    csdm.x[0] = cp.LinearDimension(
        count=10, increment="3.4 Hz", coordinates_offset="-0.05 kHz", complex_fft=True
    )
    res = {"count": 10, "spectral_width": 34, "reference_offset": -50}
    assertion(csdm, res, "even complex_fft coordinates_offset")

    # even
    csdm.x[0] = cp.LinearDimension(count=10, increment="1 Hz")
    res = {"count": 10, "spectral_width": 10, "reference_offset": 5}
    assertion(csdm, res, "even")

    csdm.x[0] = cp.LinearDimension(count=10, increment="3.4 Hz")
    res = {"count": 10, "spectral_width": 34, "reference_offset": 17}
    assertion(csdm, res, "even")

    # even # coordinates_offset
    csdm.x[0] = cp.LinearDimension(
        count=10, increment="1 Hz", coordinates_offset="5 Hz", label="1H"
    )
    res = {"count": 10, "spectral_width": 10, "reference_offset": 10, "label": "1H"}
    assertion(csdm, res, "even coordinates_offset")

    csdm.x[0] = cp.LinearDimension(
        count=10, increment="3.4 Hz", coordinates_offset="-0.05 kHz"
    )
    res = {"count": 10, "spectral_width": 34, "reference_offset": -33}
    assertion(csdm, res, "even coordinates_offset")


def test_get_spectral_dimensions_odd_count_positive_increment():
    # odd # complex_fft
    csdm = cp.as_csdm(np.arange(15))
    csdm.x[0] = cp.LinearDimension(count=15, increment="1 Hz", complex_fft=True)
    res = {"count": 15, "spectral_width": 15, "reference_offset": 0}
    assertion(csdm, res, "odd complex_fft")

    csdm.x[0] = cp.LinearDimension(count=15, increment="3.4 Hz", complex_fft=True)
    res = {"count": 15, "spectral_width": 51, "reference_offset": 0}
    assertion(csdm, res, "odd complex_fft")

    # odd # complex_fft # coordinates_offset
    csdm.x[0] = cp.LinearDimension(
        count=15, increment="1 Hz", complex_fft=True, coordinates_offset="12 Hz"
    )
    res = {"count": 15, "spectral_width": 15, "reference_offset": 12}
    assertion(csdm, res, "odd complex_fft coordinates_offset")

    csdm.x[0] = cp.LinearDimension(
        count=15, increment="3.4 Hz", coordinates_offset="12 Hz", complex_fft=True
    )
    res = {"count": 15, "spectral_width": 51, "reference_offset": 12}
    assertion(csdm, res, "odd complex_fft coordinates_offset")

    # odd
    csdm = cp.as_csdm(np.arange(15))
    csdm.x[0] = cp.LinearDimension(count=15, increment="1 Hz")
    res = {"count": 15, "spectral_width": 15, "reference_offset": 7}
    assertion(csdm, res, "odd complex_fft")

    csdm.x[0] = cp.LinearDimension(count=15, increment="3.4 Hz")
    res = {"count": 15, "spectral_width": 51, "reference_offset": 23.8}
    assertion(csdm, res, "odd complex_fft")

    # odd # coordinates_offset
    csdm.x[0] = cp.LinearDimension(
        count=15, increment="1 Hz", coordinates_offset="12 Hz"
    )
    res = {"count": 15, "spectral_width": 15, "reference_offset": 19}
    assertion(csdm, res, "odd coordinates_offset")

    csdm.x[0] = cp.LinearDimension(
        count=15, increment="3.4 Hz", coordinates_offset="12 Hz"
    )
    res = {"count": 15, "spectral_width": 51, "reference_offset": 35.8}
    assertion(csdm, res, "odd coordinates_offset")


def test_get_spectral_dimensions_even_count_negative_increment():
    # even
    csdm = cp.as_csdm(np.arange(10))
    csdm.x[0] = cp.LinearDimension(count=10, increment="-1 Hz", label="1H")
    res = {"count": 10, "spectral_width": 10, "reference_offset": -4, "label": "1H"}
    assertion(csdm, res, "even")

    csdm.x[0] = cp.LinearDimension(count=10, increment="-3.4 Hz", label="1H")
    res = {"count": 10, "spectral_width": 34, "reference_offset": -13.6, "label": "1H"}
    assertion(csdm, res, "even")

    # even # coordinates_offset
    csdm.x[0] = cp.LinearDimension(
        count=10, increment="-1 Hz", coordinates_offset="-10 Hz", label="1H"
    )
    res = {"count": 10, "spectral_width": 10, "reference_offset": -14, "label": "1H"}
    assertion(csdm, res, "even coordinates_offset")

    csdm.x[0] = cp.LinearDimension(
        count=10, increment="-3.4 Hz", coordinates_offset="-10 Hz", label="1H"
    )
    res = {"count": 10, "spectral_width": 34, "reference_offset": -23.6, "label": "1H"}
    assertion(csdm, res, "even coordinates_offset")


def test_get_spectral_dimensions_odd_count_negative_increment():
    # odd
    csdm = cp.as_csdm(np.arange(15))
    csdm.x[0] = cp.LinearDimension(count=15, increment="-1 Hz", label="1H")
    res = {"count": 15, "spectral_width": 15, "reference_offset": -7, "label": "1H"}
    assertion(csdm, res, "odd")

    csdm.x[0] = cp.LinearDimension(count=15, increment="-3.4 Hz", label="1H")
    res = {"count": 15, "spectral_width": 51, "reference_offset": -23.8, "label": "1H"}
    assertion(csdm, res, "odd")

    # even # coordinates_offset
    csdm.x[0] = cp.LinearDimension(
        count=15, increment="-1 Hz", coordinates_offset="-10 Hz", label="1H"
    )
    res = {"count": 15, "spectral_width": 15, "reference_offset": -17, "label": "1H"}
    assertion(csdm, res, "odd coordinates_offset")

    csdm.x[0] = cp.LinearDimension(
        count=15, increment="-3.4 Hz", coordinates_offset="-10 Hz", label="1H"
    )
    res = {"count": 15, "spectral_width": 51, "reference_offset": -33.8, "label": "1H"}
    assertion(csdm, res, "odd coordinates_offset")


def test_generic():
    # 7
    csdm = cp.as_csdm(np.arange(10))
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
    assertion(csdm, res, "7")
