# -*- coding: utf-8 -*-
import csdmpy as cp
import numpy as np
import pytest
from mrsimulator import signal_processing as sp

__author__ = "Maxwell C. Venetos"
__email__ = "maxvenetos@gmail.com"


def setup_read_write(operation, py_dict, fn):
    # class to dict with units
    dict_ = operation.json()
    assert dict_ == py_dict

    # read from dictionary
    b = fn.parse_dict_with_units(dict_)
    assert operation == b

    processor = sp.SignalProcessor.parse_dict_with_units({"operations": [dict_]})
    assert operation == processor.operations[0]


def test_01():
    post_sim = sp.SignalProcessor()
    operations = [
        sp.IFFT(),
        sp.apodization.Gaussian(FWHM="12 K", dim_index=0, dv_index=0),
        sp.FFT(),
    ]

    post_sim.operations = operations

    with pytest.raises(ValueError, match="The data must be a CSDM object."):
        post_sim.apply_operations([])

    data = cp.as_csdm(np.arange(20))
    data.x[0] = cp.LinearDimension(count=20, increment="10 K")
    post_sim.apply_operations(data)

    # to dict with units
    dict_ = post_sim.json()
    assert dict_ == {
        "operations": [
            {"dim_index": 0, "function": "IFFT"},
            {
                "function": "apodization",
                "type": "Gaussian",
                "FWHM": "12.0 K",
                "dim_index": 0,
                "dv_index": 0,
            },
            {"dim_index": 0, "function": "FFT"},
        ],
    }

    # parse dict with units
    post_sim_2 = sp.SignalProcessor.parse_dict_with_units(dict_)

    assert post_sim.operations == post_sim_2.operations


def generate_data():
    dv1 = cp.as_dependent_variable(np.random.rand(20))
    dv2 = cp.as_dependent_variable(np.random.rand(20))
    dv3 = cp.as_dependent_variable(np.random.rand(20))
    dim = cp.as_dimension(np.arange(20))
    return cp.CSDM(dependent_variables=[dv1, dv2, dv3], dimensions=[dim])


def test_scale():
    data_in = generate_data()
    PS_0 = [sp.Scale(factor=10)]
    operator = sp.SignalProcessor(operations=PS_0)
    data_out = operator.apply_operations(data=data_in.copy())
    _, y0, y1, y2 = data_in.to_list()
    _, y0_, y1_, y2_ = data_out.to_list()

    for in_, out_ in zip([y0, y1, y2], [y0_, y1_, y2_]):
        assert np.allclose(out_.max() / in_.max(), 10), "Scaling failed"


def test_linear():
    data_in = generate_data()
    PS_0 = [sp.Linear(amplitude=4.1, offset=10)]
    operator = sp.SignalProcessor(operations=PS_0)
    data_out = operator.apply_operations(data=data_in.copy())
    _, y0, y1, y2 = data_in.to_list()
    _, y0_, y1_, y2_ = data_out.to_list()

    for in_, out_ in zip([y0, y1, y2], [y0_, y1_, y2_]):
        assert np.allclose((out_.max() - 10) / in_.max(), 4.1), "Offset failed"
