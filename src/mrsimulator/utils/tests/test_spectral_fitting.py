# -*- coding: utf-8 -*-
import csdmpy as cp
import mrsimulator.signal_processing as sp
import mrsimulator.signal_processing.apodization as apo
import mrsimulator.utils.spectral_fitting as sf
import numpy as np
import pytest
from lmfit import Parameters
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.methods import BlochDecayCentralTransitionSpectrum


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


def test_get_and_set_simulator_object_value():
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
    sf._set_simulator_object_value(sim, string, -32)
    assert sf._get_simulator_object_value(sim, string) == -32

    string = "sys_0_site_0_shielding_symmetric_eta"
    assert sf._get_simulator_object_value(sim, string) == 0.1
    sf._set_simulator_object_value(sim, string, 0.4)
    assert sf._get_simulator_object_value(sim, string) == 0.4

    string = "sys_0_site_0_shielding_symmetric_zeta"
    assert sf._get_simulator_object_value(sim, string) == -120.0
    sf._set_simulator_object_value(sim, string, 134.2)
    assert sf._get_simulator_object_value(sim, string) == 134.2

    string = "sys_0_site_0_quadrupolar_eta"
    assert sf._get_simulator_object_value(sim, string) == 0.31
    sf._set_simulator_object_value(sim, string, 0.98)
    assert sf._get_simulator_object_value(sim, string) == 0.98

    string = "sys_0_site_0_quadrupolar_Cq"
    assert sf._get_simulator_object_value(sim, string) == 1e5
    sf._set_simulator_object_value(sim, string, 2.3e6)
    assert sf._get_simulator_object_value(sim, string) == 2.3e6

    string = "sys_0_site_0_quadrupolar_alpha"
    assert sf._get_simulator_object_value(sim, string) is None
    sf._set_simulator_object_value(sim, string, 1.32)
    assert sf._get_simulator_object_value(sim, string) == 1.32

    string = "sys_0_site_0_quadrupolar_beta"
    assert sf._get_simulator_object_value(sim, string) == 5.12
    sf._set_simulator_object_value(sim, string, None)
    assert sf._get_simulator_object_value(sim, string) is None

    string = "sys_0_site_0_quadrupolar_gamma"
    assert sf._get_simulator_object_value(sim, string) is None
    sf._set_simulator_object_value(sim, string, 1.2)
    assert sf._get_simulator_object_value(sim, string) == 1.2

    string = "sys_0_abundance"
    assert sf._get_simulator_object_value(sim, string) == 0.123
    sf._set_simulator_object_value(sim, string, 45)
    assert sf._get_simulator_object_value(sim, string) == 45


def test_03():
    # post_sim_LMFIT_params
    op_list = [
        sp.IFFT(dim_index=0),
        apo.Exponential(FWHM=100),
        apo.Gaussian(FWHM=200),
        sp.FFT(dim_index=0),
        sp.Scale(factor=10),
    ]
    post_sim = sp.SignalProcessor(operations=op_list)

    params = sf._post_sim_LMFIT_params(post_sim)

    val = params.valuesdict()

    assert val["operation_1_Exponential_FWHM"] == 100
    assert val["operation_2_Gaussian_FWHM"] == 200
    assert val["operation_4_Scale_factor"] == 10


def test_04():
    # update_post_sim_from_LMFIT_params
    op_list = [
        sp.IFFT(dim_index=0),
        apo.Gaussian(FWHM=130, dim_index=0, dv_index=0),
        apo.Exponential(FWHM=100, dim_index=1, dv_index=0),
        sp.FFT(dim_index=0),
        sp.Scale(factor=10),
    ]
    post_sim = sp.SignalProcessor(operations=op_list)

    params = sf._post_sim_LMFIT_params(post_sim)

    params["operation_1_Gaussian_FWHM"].value = 30
    params["operation_2_Exponential_FWHM"].value = 10
    params["operation_4_Scale_factor"].value = 1

    sf._update_post_sim_from_LMFIT_params(params, post_sim)

    assert post_sim.operations[1].FWHM == 30
    assert post_sim.operations[2].FWHM == 10
    assert post_sim.operations[4].factor == 1


def test_5():
    # make_LMFIT_params
    sim = Simulator()

    H = {
        "isotope": "1H",
        "isotropic_chemical_shift": "10 ppm",
        "shielding_symmetric": {"zeta": "5 ppm", "eta": 0.1},
    }
    C = {
        "isotope": "13C",
        "isotropic_chemical_shift": "-10 ppm",
        "shielding_symmetric": {"zeta": "15 ppm", "eta": 0.2},
    }
    spin_system1 = {"sites": [H], "abundance": "100%"}
    system_object1 = SpinSystem.parse_dict_with_units(spin_system1)
    spin_system2 = {"sites": [C], "abundance": "60%"}
    system_object2 = SpinSystem.parse_dict_with_units(spin_system2)
    sim.spin_systems += [system_object1, system_object2]

    op_list = [
        sp.IFFT(dim_index=0),
        apo.Exponential(FWHM=100, dim_index=0, dv_index=0),
        sp.FFT(dim_index=0),
        sp.Scale(factor=10),
    ]
    post_sim = sp.SignalProcessor(data=None, operations=op_list)

    e = "Expecting a `Simulator` object, found"
    with pytest.raises(ValueError, match=f".*{e}.*"):
        sf.make_LMFIT_params(12, 21)

    e = "Expecting a `SignalProcessor` object, found"
    with pytest.raises(ValueError, match=f".*{e}.*"):
        sf.make_LMFIT_params(sim, 21)

    params = sf.make_LMFIT_params(sim, post_sim)
    valuesdict = {
        "sys_0_site_0_isotropic_chemical_shift": 10,
        "sys_0_site_0_shielding_symmetric_zeta": 5,
        "sys_0_site_0_shielding_symmetric_eta": 0.1,
        "sys_0_abundance": 62.5,
        "sys_1_site_0_isotropic_chemical_shift": -10,
        "sys_1_site_0_shielding_symmetric_zeta": 15,
        "sys_1_site_0_shielding_symmetric_eta": 0.2,
        "sys_1_abundance": 37.5,
        "operation_1_Exponential_FWHM": 100,
        "operation_3_Scale_factor": 10,
    }
    assert params.valuesdict() == valuesdict, "Parameter creation failed"

    params = sf.make_LMFIT_params(sim)
    valuesdict = {
        "sys_0_site_0_isotropic_chemical_shift": 10,
        "sys_0_site_0_shielding_symmetric_zeta": 5,
        "sys_0_site_0_shielding_symmetric_eta": 0.1,
        "sys_0_abundance": 62.5,
        "sys_1_site_0_isotropic_chemical_shift": -10,
        "sys_1_site_0_shielding_symmetric_zeta": 15,
        "sys_1_site_0_shielding_symmetric_eta": 0.2,
        "sys_1_abundance": 37.5,
    }
    assert params.valuesdict() == valuesdict, "Parameter creation failed"

    params = sf.make_LMFIT_parameters(sim)
    valuesdict = {
        "sys_0_site_0_isotropic_chemical_shift": 10,
        "sys_0_site_0_shielding_symmetric_zeta": 5,
        "sys_0_site_0_shielding_symmetric_eta": 0.1,
        "sys_0_abundance": 62.5,
        "sys_1_site_0_isotropic_chemical_shift": -10,
        "sys_1_site_0_shielding_symmetric_zeta": 15,
        "sys_1_site_0_shielding_symmetric_eta": 0.2,
        "sys_1_abundance": 37.5,
    }
    assert params.valuesdict() == valuesdict, "Parameter creation failed"


def test_6():
    # LMFIT_min_function
    e = "Expecting a `Parameters` object, found"
    with pytest.raises(ValueError, match=f".*{e}.*"):
        _ = sf.LMFIT_min_function([], [], [])

    params = Parameters()

    e = "Expecting a `SignalProcessor` object, found"
    with pytest.raises(ValueError, match=f".*{e}.*"):
        _ = sf.LMFIT_min_function(params, [], [])

    processor = sp.SignalProcessor()
    e = "Expecting a `Simulator` object, found"
    with pytest.raises(ValueError, match=f".*{e}.*"):
        _ = sf.LMFIT_min_function(params, [], processor)

    site = Site(isotope="23Na")
    sys = SpinSystem(sites=[site], abundance=100)
    sim = Simulator()
    sim.spin_systems.append(sys)
    sim.methods = [BlochDecayCentralTransitionSpectrum(channels=["23Na"])]
    sim.methods[0].experiment = cp.as_csdm(np.zeros(1024))

    processor = sp.SignalProcessor(
        operations=[
            sp.IFFT(dim_index=0),
            apo.Gaussian(FWHM="0.2 kHz", dim_index=0),
            sp.FFT(dim_index=0),
        ]
    )

    sim.run()
    data = processor.apply_operations(sim.methods[0].simulation)

    params = sf.make_LMFIT_params(sim, processor)
    a = sf.LMFIT_min_function(params, sim, processor)
    np.testing.assert_almost_equal(-a.sum(), data.sum().real, decimal=8)
