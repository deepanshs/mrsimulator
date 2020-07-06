# -*- coding: utf-8 -*-
import mrsimulator.signal_processing as sp
import mrsimulator.signal_processing.apodization as apo


def test_01():
    post_sim = sp.SignalProcessor()
    operations = [
        sp.IFFT(),
        apo.Gaussian(sigma=12, dim_indx=0, dep_var_indx=0),
        sp.FFT(),
    ]

    post_sim.operations = operations

    # to dict with units
    dict_ = post_sim.to_dict_with_units()
    assert dict_ == {
        "operations": [
            {"dim_indx": 0, "function": "IFFT"},
            {
                "function": "apodization",
                "type": "Gaussian",
                "sigma": "12.0 Hz",
                "dim_indx": 0,
                "dep_var_indx": 0,
            },
            {"dim_indx": 0, "function": "FFT"},
        ],
    }

    # parse dict with units
    post_sim_2 = sp.SignalProcessor.parse_dict_with_units(dict_)

    assert post_sim == post_sim_2
