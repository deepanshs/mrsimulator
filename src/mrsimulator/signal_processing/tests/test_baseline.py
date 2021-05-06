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
    x_c = data_in.dimensions[0].coordinates.to("s").value

    # zeroth order
    PS_0 = [sp.baseline.Polynomial(x0="10")]
    data_out = sp.SignalProcessor(operations=PS_0).apply_operations(data=data_in.copy())

    _, y0, y1, y2 = data_in.to_list()
    _, y0_, y1_, y2_ = data_out.to_list()

    for in_, out_ in zip([y0, y1, y2], [y0_, y1_, y2_]):
        assert np.allclose(out_.max() - in_.max(), 10), "Offset failed"

    py_dict = {
        "function": "baseline",
        "type": "Polynomial",
        "x0": "10.0",
        "x1": "0",
        "x2": "0",
        "x3": "0",
        "dim_index": 0,
    }
    setup_read_write(PS_0[0], py_dict, sp.baseline.Polynomial)

    # first order
    PS_0 = [sp.baseline.Polynomial(x0="30", x1="1 Hz")]
    data_out = sp.SignalProcessor(operations=PS_0).apply_operations(data=data_in.copy())

    _, y0, y1, y2 = data_in.to_list()
    _, y0_, y1_, y2_ = data_out.to_list()

    for in_, out_ in zip([y0, y1, y2], [y0_, y1_, y2_]):
        assert np.allclose(out_, in_ + 30 + x_c), "Polynomial 1st order failed"

    py_dict = {
        "function": "baseline",
        "type": "Polynomial",
        "x0": 30.0,
        "x1": "1.0 Hz",
        "x2": 0,
        "x3": 0,
        "dim_index": 0,
    }
    setup_read_write(PS_0[0], py_dict, sp.baseline.Polynomial)

    # second order
    PS_0 = [sp.baseline.Polynomial(x0="1", x2="1 Hz^2")]
    data_out = sp.SignalProcessor(operations=PS_0).apply_operations(data=data_in.copy())

    _, y0, y1, y2 = data_in.to_list()
    _, y0_, y1_, y2_ = data_out.to_list()

    for in_, out_ in zip([y0, y1, y2], [y0_, y1_, y2_]):
        assert np.allclose(out_, in_ + 1 + x_c ** 2), "Polynomial 2nd order failed"

    py_dict = {
        "function": "baseline",
        "type": "Polynomial",
        "x0": 1.0,
        "x1": 0,
        "x2": "1.0 Hz2",
        "x3": 0,
        "dim_index": 0,
    }
    setup_read_write(PS_0[0], py_dict, sp.baseline.Polynomial)

    # third order
    PS_0 = [sp.baseline.Polynomial(x0="10", x3="2 Hz^3", x1="13.1 Hz")]
    data_out = sp.SignalProcessor(operations=PS_0).apply_operations(data=data_in.copy())

    _, y0, y1, y2 = data_in.to_list()
    _, y0_, y1_, y2_ = data_out.to_list()

    for in_, out_ in zip([y0, y1, y2], [y0_, y1_, y2_]):
        assert np.allclose(
            out_, in_ + 10 + 13.1 * x_c + 2 * x_c ** 3
        ), "Polynomial 3rd order failed"

    # py_dict = {
    #     "function": "baseline",
    #     "type": "Polynomial",
    #     "x0": 10.0,
    #     "x1": "1.0 Hz",
    #     "x2": 0,
    #     "x3": "2.0 Hz3",
    #     "dim_index": 0,
    # }
    # setup_read_write(PS_0[0], py_dict, sp.baseline.Polynomial)
