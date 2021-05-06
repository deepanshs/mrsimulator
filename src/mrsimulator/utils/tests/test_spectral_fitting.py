# -*- coding: utf-8 -*-
import csdmpy as cp
import mrsimulator.utils.spectral_fitting as sf
import numpy as np
import pytest
from mrsimulator import signal_processing as sp
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.methods import BlochDecayCTSpectrum

# from lmfit import Parameters


def test_str_encode():
    str1 = "spin_systems[0].sites[0].isotropic_chemical_shift"
    str1_encoded = "sys_0_site_0_isotropic_chemical_shift"

    str2 = "spin_systems[0].sites[0].quadrupolar.Cq"
    str2_encoded = "sys_0_site_0_quadrupolar_Cq"

    assert sf._str_encode(str1) == str1_encoded
    assert sf._str_encode(str2) == str2_encoded


def test_str_decode():
    str1 = ["spin_systems", "0", "sites", "0", "isotropic_chemical_shift"]
    str1_encoded = "sys_0_site_0_isotropic_chemical_shift"

    str2 = ["spin_systems", "0", "sites", "0", "quadrupolar", "Cq"]
    str2_encoded = "sys_0_site_0_quadrupolar_Cq"

    assert sf._str_decode(str1_encoded) == str1
    assert sf._str_decode(str2_encoded) == str2


def test_get_and_update_sim_from_LMFIT_params():
    site = Site(
        isotope="23Na",
        isotropic_chemical_shift=32,
        shielding_symmetric={"zeta": -120, "eta": 0.1},
        quadrupolar={"Cq": 1e5, "eta": 0.31, "beta": 5.12},
    )
    sys = SpinSystem(sites=[site], abundance=0.123)
    sim = Simulator()
    sim.spin_systems.append(sys)

    string = "sys_0_site_0_isotropic_chemical_shift"
    assert sf._get_simulator_object_value(sim, string) == 32
    sf._update_sim_from_LMFIT_params(sim, string, -32)
    assert sf._get_simulator_object_value(sim, string) == -32

    string = "sys_0_site_0_shielding_symmetric_eta"
    assert sf._get_simulator_object_value(sim, string) == 0.1
    sf._update_sim_from_LMFIT_params(sim, string, 0.4)
    assert sf._get_simulator_object_value(sim, string) == 0.4

    string = "sys_0_site_0_shielding_symmetric_zeta"
    assert sf._get_simulator_object_value(sim, string) == -120.0
    sf._update_sim_from_LMFIT_params(sim, string, 134.2)
    assert sf._get_simulator_object_value(sim, string) == 134.2

    string = "sys_0_site_0_quadrupolar_eta"
    assert sf._get_simulator_object_value(sim, string) == 0.31
    sf._update_sim_from_LMFIT_params(sim, string, 0.98)
    assert sf._get_simulator_object_value(sim, string) == 0.98

    string = "sys_0_site_0_quadrupolar_Cq"
    assert sf._get_simulator_object_value(sim, string) == 1e5
    sf._update_sim_from_LMFIT_params(sim, string, 2.3e6)
    assert sf._get_simulator_object_value(sim, string) == 2.3e6

    string = "sys_0_site_0_quadrupolar_alpha"
    assert sf._get_simulator_object_value(sim, string) is None
    sf._update_sim_from_LMFIT_params(sim, string, 1.32)
    assert sf._get_simulator_object_value(sim, string) == 1.32

    string = "sys_0_site_0_quadrupolar_beta"
    assert sf._get_simulator_object_value(sim, string) == 5.12
    sf._update_sim_from_LMFIT_params(sim, string, None)
    assert sf._get_simulator_object_value(sim, string) is None

    string = "sys_0_site_0_quadrupolar_gamma"
    assert sf._get_simulator_object_value(sim, string) is None
    sf._update_sim_from_LMFIT_params(sim, string, 1.2)
    assert sf._get_simulator_object_value(sim, string) == 1.2

    string = "sys_0_abundance"
    assert sf._get_simulator_object_value(sim, string) == 0.123
    sf._update_sim_from_LMFIT_params(sim, string, 45)
    assert sf._get_simulator_object_value(sim, string) == 45


def test_03():
    # post_sim_LMFIT_params
    op_list = [
        sp.IFFT(dim_index=0),
        sp.apodization.Exponential(FWHM=100),
        sp.apodization.Gaussian(FWHM=200),
        sp.FFT(dim_index=0),
        sp.Scale(factor=10),
        sp.baseline.ConstantOffset(offset=43.1),
        sp.Linear(amplitude=32.9, offset=13.4),
    ]
    post_sim = sp.SignalProcessor(operations=op_list)

    params = sf.make_signal_processor_params(post_sim)

    val = params.valuesdict()

    assert val["SP_0_operation_1_Exponential_FWHM"] == 100
    assert val["SP_0_operation_2_Gaussian_FWHM"] == 200
    assert val["SP_0_operation_4_Scale_factor"] == 10
    assert val["SP_0_operation_5_ConstantOffset_offset"] == 43.1
    assert val["SP_0_operation_6_Linear_amplitude"] == 32.9
    assert val["SP_0_operation_6_Linear_offset"] == 13.4


def test_04():
    # update_post_sim_from_LMFIT_params
    op_list1 = [sp.Scale(factor=20)]
    op_list2 = [
        sp.IFFT(dim_index=0),
        sp.apodization.Gaussian(FWHM=130, dim_index=0, dv_index=0),
        sp.apodization.Exponential(FWHM=100, dim_index=1, dv_index=0),
        sp.FFT(dim_index=0),
        sp.Scale(factor=10),
    ]
    post_sim = [sp.SignalProcessor(operations=lst) for lst in [op_list1, op_list2]]

    params = sf.make_signal_processor_params(post_sim)

    params["SP_0_operation_0_Scale_factor"].value = 300
    params["SP_1_operation_1_Gaussian_FWHM"].value = 30
    params["SP_1_operation_2_Exponential_FWHM"].value = 10
    params["SP_1_operation_4_Scale_factor"].value = 1

    sf._update_post_sim_from_LMFIT_params(params, post_sim)

    assert post_sim[0].operations[0].factor == 300
    assert post_sim[1].operations[1].FWHM == 30
    assert post_sim[1].operations[2].FWHM == 10
    assert post_sim[1].operations[4].factor == 1


def test_raise_messages():
    e = "Expecting a `Simulator` object, found"
    with pytest.raises(ValueError, match=f".*{e}.*"):
        sf.make_LMFIT_params(12, 21)

    e = "Expecting a `SignalProcessor` object, found"
    with pytest.raises(ValueError, match=f".*{e}.*"):
        sf.make_LMFIT_params(Simulator(spin_systems=[SpinSystem()]), 21)


H = {
    "isotope": "1H",
    "isotropic_chemical_shift": "10 ppm",
    "shielding_symmetric": {
        "zeta": "5 ppm",
        "eta": 0.1,
        "alpha": "3.12 rad",
        "gamma": "0.341 rad",
    },
}
C = {
    "isotope": "13C",
    "isotropic_chemical_shift": "-10 ppm",
    "shielding_symmetric": {"zeta": "15 ppm", "eta": 0.2, "beta": "4.12 rad"},
}
CH = {
    "site_index": [0, 1],
    "isotropic_j": "10 Hz",
    "j_symmetric": {"zeta": "60 Hz", "eta": 0.4},
}
op_list = [
    sp.IFFT(dim_index=0),
    sp.apodization.Exponential(FWHM=100, dim_index=0, dv_index=0),
    sp.FFT(dim_index=0),
    sp.Scale(factor=10),
]


def test_5_multi_spin_systems():
    sim = Simulator()
    spin_system1 = {"sites": [H], "abundance": "100%"}
    system_object1 = SpinSystem.parse_dict_with_units(spin_system1)
    spin_system2 = {"sites": [C], "abundance": "60%"}
    system_object2 = SpinSystem.parse_dict_with_units(spin_system2)
    sim.spin_systems += [system_object1, system_object2]
    post_sim = sp.SignalProcessor(operations=op_list)

    params = sf.make_LMFIT_params(sim, post_sim)
    valuesdict_system = {
        "sys_0_site_0_isotropic_chemical_shift": 10,
        "sys_0_site_0_shielding_symmetric_zeta": 5,
        "sys_0_site_0_shielding_symmetric_eta": 0.1,
        "sys_0_site_0_shielding_symmetric_alpha": 3.12,
        "sys_0_site_0_shielding_symmetric_gamma": 0.341,
        "sys_0_abundance": 62.5,
        "sys_1_site_0_isotropic_chemical_shift": -10,
        "sys_1_site_0_shielding_symmetric_zeta": 15,
        "sys_1_site_0_shielding_symmetric_eta": 0.2,
        "sys_1_site_0_shielding_symmetric_beta": 4.12,
        "sys_1_abundance": 37.5,
    }
    valuedict_proc = {
        "SP_0_operation_1_Exponential_FWHM": 100,
        "SP_0_operation_3_Scale_factor": 10,
    }
    assert params.valuesdict() == {
        **valuesdict_system,
        **valuedict_proc,
    }, "Parameter creation failed"

    params = sf.make_LMFIT_params(sim)
    assert params.valuesdict() == valuesdict_system, "Parameter creation failed"

    # alias
    params = sf.make_LMFIT_parameters(sim)
    assert params.valuesdict() == valuesdict_system, "Parameter creation failed"


def test_6_coupled():
    sim = Simulator()
    spin_system = {"sites": [H, C], "couplings": [CH], "abundance": "100%"}
    system_object = SpinSystem.parse_dict_with_units(spin_system)
    sim.spin_systems += [system_object]
    post_sim = sp.SignalProcessor(operations=op_list)

    params = sf.make_LMFIT_params(sim, post_sim)
    valuesdict_system = {
        "sys_0_site_0_isotropic_chemical_shift": 10,
        "sys_0_site_0_shielding_symmetric_zeta": 5,
        "sys_0_site_0_shielding_symmetric_eta": 0.1,
        "sys_0_site_0_shielding_symmetric_alpha": 3.12,
        "sys_0_site_0_shielding_symmetric_gamma": 0.341,
        "sys_0_site_1_isotropic_chemical_shift": -10,
        "sys_0_site_1_shielding_symmetric_zeta": 15,
        "sys_0_site_1_shielding_symmetric_eta": 0.2,
        "sys_0_site_1_shielding_symmetric_beta": 4.12,
        "sys_0_coupling_0_isotropic_j": 10,
        "sys_0_coupling_0_j_symmetric_zeta": 60,
        "sys_0_coupling_0_j_symmetric_eta": 0.4,
        "sys_0_abundance": 100,
    }
    valuedict_proc = {
        "SP_0_operation_1_Exponential_FWHM": 100,
        "SP_0_operation_3_Scale_factor": 10,
    }
    assert params.valuesdict() == {
        **valuesdict_system,
        **valuedict_proc,
    }, "Parameter creation failed"

    params = sf.make_LMFIT_params(sim)
    assert params.valuesdict() == valuesdict_system, "Parameter creation failed"

    # alias
    params = sf.make_LMFIT_parameters(sim)
    assert params.valuesdict() == valuesdict_system, "Parameter creation failed"


def test_7():
    site = Site(isotope="23Na")
    sys = SpinSystem(sites=[site], abundance=50)
    sim = Simulator()
    sim.spin_systems = [sys, sys]
    sim.methods = [BlochDecayCTSpectrum(channels=["23Na"])]
    sim.methods[0].experiment = cp.as_csdm(np.zeros(1024))

    processor = sp.SignalProcessor(
        operations=[
            sp.IFFT(dim_index=0),
            sp.apodization.Gaussian(FWHM="0.2 kHz", dim_index=0),
            sp.FFT(dim_index=0),
        ]
    )

    def test_array():
        sim.run()
        data = processor.apply_operations(sim.methods[0].simulation)

        data_sum = 0
        for dv in data.y:
            data_sum += dv.components[0]

        params = sf.make_LMFIT_params(sim, processor)
        a = sf.LMFIT_min_function(params, sim, processor)
        np.testing.assert_almost_equal(-a, data_sum, decimal=8)

    test_array()

    sim.config.decompose_spectrum = "spin_system"
    test_array()

    # data = processor.apply_operations(sim.methods[0].simulation)

    # params = sf.make_LMFIT_params(sim, processor)
    # a = sf.LMFIT_min_function(params, sim, processor)
    # np.testing.assert_almost_equal(-a.sum(), data.sum().real, decimal=8)
