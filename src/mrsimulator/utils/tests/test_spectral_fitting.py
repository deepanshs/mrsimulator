# -*- coding: utf-8 -*-
import csdmpy as cp
import mrsimulator.utils.spectral_fitting as sf
import numpy as np
import pytest
from mrsimulator import signal_processing as sp
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.method.lib import BlochDecayCTSpectrum
from mrsimulator.method.lib import BlochDecaySpectrum
from mrsimulator.method.lib import SSB2D
from mrsimulator.method.lib import ThreeQ_VAS


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


def setup_simulator():
    site = Site(
        isotope="23Na",
        isotropic_chemical_shift=32,
        shielding_symmetric={"zeta": -120, "eta": 0.1},
        quadrupolar={"Cq": 1e5, "eta": 0.31, "beta": 5.12},
    )
    sys = SpinSystem(sites=[site], abundance=0.123)
    sim = Simulator()
    sim.spin_systems.append(sys)
    sim.methods.append(BlochDecayCTSpectrum(channels=["2H"], rotor_frequency=1e3))
    sim.methods.append(BlochDecaySpectrum(channels=["2H"], rotor_frequency=12.5e3))
    sim.methods.append(ThreeQ_VAS(channels=["27Al"]))
    sim.methods.append(SSB2D(channels=["17O"], rotor_frequency=35500))
    return sim


def test_make_simulator_params():
    sim = setup_simulator()

    def check(val):
        assert val["sys_0_site_0_isotropic_chemical_shift"] == 32
        assert val["sys_0_site_0_shielding_symmetric_zeta"] == -120
        assert val["sys_0_site_0_shielding_symmetric_eta"] == 0.1
        assert val["sys_0_site_0_quadrupolar_Cq"] == 1e5
        assert val["sys_0_site_0_quadrupolar_eta"] == 0.31
        assert val["sys_0_site_0_quadrupolar_beta"] == 5.12
        assert val["sys_0_abundance"] == 100

    params = sf.make_simulator_params(sim)
    check(params.valuesdict())

    params = sf.make_simulator_params(sim, include={"rotor_frequency"})
    val = params.valuesdict()
    check(params.valuesdict())
    assert val["mth_0_rotor_frequency"] == 1e3
    assert val["mth_1_rotor_frequency"] == 12.5e3
    assert val["mth_3_rotor_frequency"] == 35500


def test_update_simulator_from_LMFIT_params():
    sim = setup_simulator()

    def set_params(params):
        params["sys_0_site_0_isotropic_chemical_shift"].value = -110
        params["sys_0_site_0_shielding_symmetric_zeta"].value = 50
        params["sys_0_site_0_shielding_symmetric_eta"].value = 0.9
        params["sys_0_site_0_quadrupolar_Cq"].value = 5e6
        params["sys_0_site_0_quadrupolar_eta"].value = 0.1
        params["sys_0_site_0_quadrupolar_beta"].value = 2.12
        params["sys_0_abundance"].value = 20

    def check_object(sim):
        assert sim.spin_systems[0].sites[0].isotropic_chemical_shift == -110
        assert sim.spin_systems[0].sites[0].shielding_symmetric.zeta == 50
        assert sim.spin_systems[0].sites[0].shielding_symmetric.eta == 0.9
        assert sim.spin_systems[0].sites[0].quadrupolar.Cq == 5e6
        assert sim.spin_systems[0].sites[0].quadrupolar.eta == 0.1
        assert sim.spin_systems[0].sites[0].quadrupolar.beta == 2.12
        assert sim.spin_systems[0].abundance == 100

    params = sf.make_simulator_params(sim)
    set_params(params)
    sf._update_simulator_from_LMFIT_params(params, sim)
    check_object(sim)

    params = sf.make_simulator_params(sim, include={"rotor_frequency"})
    set_params(params)
    params["mth_0_rotor_frequency"].value = 1051.120
    params["mth_1_rotor_frequency"].value = 12442.12
    params["mth_3_rotor_frequency"].value = 35550
    sf._update_simulator_from_LMFIT_params(params, sim)
    check_object(sim)
    assert sim.methods[0].spectral_dimensions[0].events[0].rotor_frequency == 1051.120
    assert sim.methods[1].spectral_dimensions[0].events[0].rotor_frequency == 12442.12
    assert sim.methods[2].spectral_dimensions[0].events[0].rotor_frequency == 1e12
    assert sim.methods[2].spectral_dimensions[1].events[0].rotor_frequency == 1e12
    assert sim.methods[3].spectral_dimensions[0].events[0].rotor_frequency == 35550
    assert sim.methods[3].spectral_dimensions[1].events[0].rotor_frequency == 1e12


def setup_signal_processor():
    op_list1 = [
        sp.IFFT(dim_index=0),
        sp.apodization.Exponential(FWHM=100),
        sp.apodization.Gaussian(FWHM=200),
        sp.FFT(dim_index=0),
        sp.Scale(factor=10),
        sp.baseline.ConstantOffset(offset=43.1),
        sp.Linear(amplitude=32.9, offset=13.4),
    ]
    op_list2 = [
        sp.Scale(factor=20),
        sp.baseline.ConstantOffset(offset=-43.1),
        sp.Linear(amplitude=1.2, offset=0.4),
    ]
    return [sp.SignalProcessor(operations=op) for op in [op_list1, op_list2]]


def test_make_signal_processor_params():
    post_sim = setup_signal_processor()
    params = sf.make_signal_processor_params(post_sim)

    val = params.valuesdict()
    assert val["SP_0_operation_1_Exponential_FWHM"] == 100
    assert val["SP_0_operation_2_Gaussian_FWHM"] == 200
    assert val["SP_0_operation_4_Scale_factor"] == 10
    assert val["SP_0_operation_5_ConstantOffset_offset"] == 43.1
    assert val["SP_0_operation_6_Linear_amplitude"] == 32.9
    assert val["SP_0_operation_6_Linear_offset"] == 13.4

    assert val["SP_1_operation_0_Scale_factor"] == 20
    assert val["SP_1_operation_1_ConstantOffset_offset"] == -43.1
    assert val["SP_1_operation_2_Linear_amplitude"] == 1.2
    assert val["SP_1_operation_2_Linear_offset"] == 0.4


def test_update_processors_from_LMFIT_params():
    post_sim = setup_signal_processor()
    params = sf.make_signal_processor_params(post_sim)

    params["SP_0_operation_1_Exponential_FWHM"].value = 30
    params["SP_0_operation_2_Gaussian_FWHM"].value = 10
    params["SP_0_operation_4_Scale_factor"].value = 1
    params["SP_0_operation_5_ConstantOffset_offset"].value = 2
    params["SP_0_operation_6_Linear_amplitude"].value = 1.2
    params["SP_0_operation_6_Linear_offset"].value = -0.1

    params["SP_1_operation_0_Scale_factor"].value = -20
    params["SP_1_operation_1_ConstantOffset_offset"].value = 13.1
    params["SP_1_operation_2_Linear_amplitude"].value = 5.1
    params["SP_1_operation_2_Linear_offset"].value = -40.1

    sf._update_processors_from_LMFIT_params(params, post_sim)

    assert post_sim[0].operations[1].FWHM == 30
    assert post_sim[0].operations[2].FWHM == 10
    assert post_sim[0].operations[4].factor == 1
    assert post_sim[0].operations[5].offset == 2
    assert post_sim[0].operations[6].amplitude == 1.2
    assert post_sim[0].operations[6].offset == -0.1

    assert post_sim[1].operations[0].factor == -20
    assert post_sim[1].operations[1].offset == 13.1
    assert post_sim[1].operations[2].amplitude == 5.1
    assert post_sim[1].operations[2].offset == -40.1


def test_raise_messages():
    e = "Expecting a `Simulator` object, found"
    with pytest.raises(ValueError, match=f".*{e}.*"):
        sf.make_LMFIT_params(12, 21)

    e = "Expecting a list of `SignalProcessor` objects."
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


def compare_result(params, valuesdict_system, sim):
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
    compare_result(params, valuesdict_system, sim)


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
    compare_result(params, valuesdict_system, sim)


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

        dat = sf.add_csdm_dvs(data.real)
        fits = sf.bestfit(sim, processor)
        assert sf.add_csdm_dvs(fits[0]) == dat

        res = sf.residuals(sim, processor)
        assert res[0] == -dat

    test_array()

    sim.config.decompose_spectrum = "spin_system"
    test_array()

    # data = processor.apply_operations(sim.methods[0].simulation)

    # params = sf.make_LMFIT_params(sim, processor)
    # a = sf.LMFIT_min_function(params, sim, processor)
    # np.testing.assert_almost_equal(-a.sum(), data.sum().real, decimal=8)
