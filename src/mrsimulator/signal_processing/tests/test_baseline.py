# -*- coding: utf-8 -*-
import csdmpy as cp
import numpy as np
from mrsimulator import signal_processing as sp

from .test_signal_processing import generate_data
from .test_signal_processing import setup_read_write


def test_baseline_constant_offset():
    data_in = generate_data()
    PS_0 = [sp.baseline.ConstantOffset(offset=10)]
    operator = sp.SignalProcessor(operations=PS_0)
    data_out = operator.apply_operations(data=data_in.copy())
    _, y0, y1, y2 = data_in.to_list()
    _, y0_, y1_, y2_ = data_out.to_list()

    for in_, out_ in zip([y0, y1, y2], [y0_, y1_, y2_]):
        assert np.allclose(out_.max() - in_.max(), 10), "Offset failed"

    py_dict = {
        "function": "baseline",
        "type": "ConstantOffset",
        "offset": 10.0,
    }
    setup_read_write(PS_0[0], py_dict, sp.baseline.ConstantOffset)


def test_baseline_polynomial():
    data_in = generate_data()
    data_in.dimensions[0] *= cp.ScalarQuantity("1 ms")
    x_c = data_in.dimensions[0].coordinates.value

    # zeroth order
    PS_0 = [sp.baseline.Polynomial(polynomial_dictionary={"c0": 10})]
    data_out = sp.SignalProcessor(operations=PS_0).apply_operations(data=data_in.copy())

    _, y0, y1, y2 = data_in.to_list()
    _, y0_, y1_, y2_ = data_out.to_list()

    for in_, out_ in zip([y0, y1, y2], [y0_, y1_, y2_]):
        assert np.allclose(out_.max() - in_.max(), 10), "Offset failed"

    # first order
    PS_0 = [sp.baseline.Polynomial(polynomial_dictionary={"c0": 30, "c1": 1})]
    data_out = sp.SignalProcessor(operations=PS_0).apply_operations(data=data_in.copy())

    _, y0, y1, y2 = data_in.to_list()
    _, y0_, y1_, y2_ = data_out.to_list()

    for in_, out_ in zip([y0, y1, y2], [y0_, y1_, y2_]):
        assert np.allclose(out_, in_ + 30 + x_c), "Polynomial 1st order failed"

    # second order
    PS_0 = [sp.baseline.Polynomial(polynomial_dictionary={"c0": 1, "c2": 1})]
    data_out = sp.SignalProcessor(operations=PS_0).apply_operations(data=data_in.copy())

    _, y0, y1, y2 = data_in.to_list()
    _, y0_, y1_, y2_ = data_out.to_list()

    for in_, out_ in zip([y0, y1, y2], [y0_, y1_, y2_]):
        assert np.allclose(out_, in_ + 1 + x_c**2), "Polynomial 2nd order failed"

    # third order
    PS_0 = [
        sp.baseline.Polynomial(polynomial_dictionary={"c0": 10, "c3": 2, "c1": 13.1})
    ]
    data_out = sp.SignalProcessor(operations=PS_0).apply_operations(data=data_in.copy())

    _, y0, y1, y2 = data_in.to_list()
    _, y0_, y1_, y2_ = data_out.to_list()

    for in_, out_ in zip([y0, y1, y2], [y0_, y1_, y2_]):
        assert np.allclose(
            out_, in_ + 10 + 13.1 * x_c + 2 * x_c**3
        ), "Polynomial3rd order failed"
