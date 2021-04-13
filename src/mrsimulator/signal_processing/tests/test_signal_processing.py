# -*- coding: utf-8 -*-
import csdmpy as cp
import numpy as np
import pytest
from mrsimulator import signal_processing as sp

__author__ = "Maxwell C. Venetos"
__email__ = "maxvenetos@gmail.com"


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
