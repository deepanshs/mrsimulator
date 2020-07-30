# -*- coding: utf-8 -*-
import mrsimulator.signal_processing as sp
import mrsimulator.signal_processing.apodization as apo
import mrsimulator.utils.spectral_fitting as sf
from mrsimulator import Simulator
from mrsimulator import SpinSystem


def test_01():
    # str_to_html
    str1 = "spin_systems[0].sites[0].isotropic_chemical_shift"
    str1_encoded = "sys_0_site_0_isotropic_chemical_shift"

    str2 = "spin_systems[0].sites[0].quadrupolar.Cq"
    str2_encoded = "sys_0_site_0_quadrupolar_Cq"

    assert sf._str_encode(str1) == str1_encoded
    assert sf._str_encode(str2) == str2_encoded


def test_02():
    # html_to_str
    str1 = "spin_systems[0].sites[0].isotropic_chemical_shift"
    str1_encoded = "sys_0_site_0_isotropic_chemical_shift"

    str2 = "spin_systems[0].sites[0].quadrupolar.Cq"
    str2_encoded = "sys_0_site_0_quadrupolar_Cq"

    assert sf._str_decode(str1_encoded) == str1
    assert sf._str_decode(str2_encoded) == str2


def test_03():
    # post_sim_LMFIT_params
    op_list = [
        sp.IFFT(dim_indx=0),
        apo.Exponential(FWHM=100, dim_indx=0, dv_indx=0),
        sp.FFT(dim_indx=0),
        sp.Scale(factor=10),
    ]
    post_sim = sp.SignalProcessor(data=None, operations=op_list)

    params = sf._post_sim_LMFIT_params(post_sim)

    val = params.valuesdict()

    assert val["operation_1_Exponential_FWHM"] == 100
    assert val["operation_3_Scale"] == 10


def test_04():
    # update_post_sim_from_LMFIT_params
    op_list = [
        sp.IFFT(dim_indx=0),
        apo.Exponential(FWHM=100, dim_indx=0, dv_indx=0),
        sp.FFT(dim_indx=0),
        sp.Scale(factor=10),
    ]
    post_sim = sp.SignalProcessor(data=None, operations=op_list)

    params = sf._post_sim_LMFIT_params(post_sim)

    params["operation_1_Exponential_FWHM"].value = 10
    params["operation_3_Scale"].value = 1

    sf._update_post_sim_from_LMFIT_params(params, post_sim)

    assert post_sim.operations[1].FWHM == 10
    assert post_sim.operations[3].factor == 1


def test_5():
    # make_LMFIT_parameters
    sim = Simulator()

    H = {
        "isotope": "1H",
        "isotropic_chemical_shift": "10 ppm",
        "shielding_symmetric": {"zeta": "5 ppm", "eta": 0.1},
    }
    spin_system = {"sites": [H], "abundance": "100%"}
    system_object = SpinSystem.parse_dict_with_units(spin_system)
    sim.spin_systems += [system_object]

    op_list = [
        sp.IFFT(dim_indx=0),
        apo.Exponential(FWHM=100, dim_indx=0, dv_indx=0),
        sp.FFT(dim_indx=0),
        sp.Scale(factor=10),
    ]
    post_sim = sp.SignalProcessor(data=None, operations=op_list)

    params = sf.make_LMFIT_parameters(sim, post_sim)

    valuesdict = {
        "sys_0_site_0_isotropic_chemical_shift": 10,
        "sys_0_site_0_shielding_symmetric_zeta": 5,
        "sys_0_site_0_shielding_symmetric_eta": 0.1,
        "sys_0_abundance": 100,
        "operation_1_Exponential_FWHM": 100,
        "operation_3_Scale": 10,
    }

    assert params.valuesdict() == valuesdict, "Parameter creation failed"
