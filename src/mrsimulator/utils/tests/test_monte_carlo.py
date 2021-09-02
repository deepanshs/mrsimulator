# -*- coding: utf-8 -*-
import copy

import csdmpy as cp
import mrsimulator.utils.monte_carlo as mc
import mrsimulator.utils.spectral_fitting as sf
import numpy as np
import pytest  # noqa
from mrsimulator import signal_processing as sp
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.methods import BlochDecayCTSpectrum
from mrsimulator.methods import BlochDecaySpectrum


def test_name_abbrev():
    params = setup_params()

    abv_0 = [
        "delta_iso_0_0",
        "zeta_0_0",
        "etaCS_0_0",
        "Cq_0_0",
        "etaQ_0_0",
        "abundance_0",
        "spin_freq_0",
        "spin_freq_1",
        "expo_fwhm_0_1",
        "gauss_fwhm_0_2",
        "scale_0_4",
        "offset_0_5",
        "linear_0_6",
        "scale_1_0",
        "offset_1_1",
        "linear_1_2",
    ]

    abv_1 = mc.name_abbrev(params)

    assert abv_0 == abv_1


def test_vary_only():
    params = setup_params()

    params["sys_0_abundance"].vary = True
    params["SP_0_operation_1_Exponential_FWHM"].vary = False

    variables = copy.deepcopy(params)
    variables.pop("sys_0_abundance")
    variables.pop("SP_0_operation_1_Exponential_FWHM")

    assert mc.mrsim_emcee._vary_only(params) == variables


def test_log_prior():
    maximum = [1, np.inf]
    minimum = [0, -np.inf]
    bounds = np.array([maximum, minimum])
    theta_1 = [0.5, 10]
    theta_2 = [10, 10]

    assert mc.mrsim_emcee._log_prior(theta_1, bounds) == 0.0
    assert mc.mrsim_emcee._log_prior(theta_2, bounds) == -np.inf


def test_log_probability():
    sim = setup_simulator()
    processor = setup_signal_processor()
    params = sf.make_LMFIT_params(sim, processor, include="rotor_frequency")
    theta_in = np.zeros(len(params)) + 0.5
    theta_in[7] = 1000
    theta_in[8] = 12500
    theta_out = np.zeros(len(params)) + 5
    emcee_obj = mc.mrsim_emcee(params, sim, processor)

    log_prob = emcee_obj._log_probability(theta_in, params, sim, processor, None)
    assert log_prob == -2304.0

    log_prob = emcee_obj._log_probability(theta_out, params, sim, processor, None)
    assert log_prob == -np.inf


def test_update_sim_params():
    sim_2 = setup_second_simulator()
    params_1 = setup_params()
    sim_2 = mc.mrsim_emcee._update_sim_params(sim_2, params_1)

    assert sim_2.spin_systems[0].sites[0].isotropic_chemical_shift == 32.0
    assert sim_2.spin_systems[0].sites[0].shielding_symmetric.zeta == -120
    assert sim_2.spin_systems[0].sites[0].shielding_symmetric.eta == 0.1
    assert sim_2.spin_systems[0].sites[0].quadrupolar.Cq == 1e5
    assert sim_2.spin_systems[0].sites[0].quadrupolar.eta == 0.31
    assert sim_2.spin_systems[0].sites[0].quadrupolar.beta == 5.12


def test_update_methods():
    sim_2 = setup_second_simulator()
    params_1 = setup_params()
    sim_2 = mc.mrsim_emcee._update_methods(sim_2, params_1)

    assert sim_2.methods[0].spectral_dimensions[0].events[0].rotor_frequency == 1000.0
    assert sim_2.methods[1].spectral_dimensions[0].events[0].rotor_frequency == 12500.0


def test_update_signal_processors():
    processor_1 = setup_signal_processor()
    processor_2 = setup_second_signal_processor()
    params_1 = setup_params()
    processor_2 = mc.mrsim_emcee._update_signal_processors(processor_2, params_1)

    for p in enumerate(processor_1):
        assert processor_1[p[0]].operations == processor_2[p[0]].operations


def test_mcmc():
    sim = setup_simulator()
    processor = setup_signal_processor()
    params = sf.make_LMFIT_params(sim, processor, include="rotor_frequency")
    emcee_obj = mc.mrsim_emcee(params, sim, processor)

    result = emcee_obj.mcmc(steps=50, nwalkers=50, burn=10, thin=5)
    assert ((result["accept_frac"] < 1) & (result["accept_frac"] > 0)).sum() == 50


def setup_simulator():
    site = Site(
        isotope="23Na",
        isotropic_chemical_shift=32,
        shielding_symmetric={"zeta": -120, "eta": 0.1},
        quadrupolar={"Cq": 1e5, "eta": 0.31, "beta": 5.12},
    )
    sys = SpinSystem(sites=[site], abundance=0.123)
    sim = Simulator()
    spectral_dims = [
        {
            "count": 4096,
            "spectral_width": 20000.0,
            "reference_offset": 1220.7031,
            "origin_offset": 54237080.0,
        }
    ]
    experiment = cp.as_csdm(np.zeros(4096))
    sim.spin_systems.append(sys)
    sim.methods.append(
        BlochDecayCTSpectrum(
            channels=["2H"],
            rotor_frequency=1e3,
            spectral_dimensions=spectral_dims,
            experiment=experiment,
        )
    )
    sim.methods.append(
        BlochDecaySpectrum(
            channels=["2H"],
            rotor_frequency=12.5e3,
            spectral_dimensions=spectral_dims,
            experiment=experiment,
        )
    )
    return sim


def setup_second_simulator():
    site = Site(
        isotope="27Al",
        isotropic_chemical_shift=40,
        shielding_symmetric={"zeta": -100, "eta": 0.45},
        quadrupolar={"Cq": 1e6, "eta": 0.31, "beta": 6.37},
    )
    sys = SpinSystem(sites=[site], abundance=0.455)
    sim = Simulator()
    spectral_dims = [
        {
            "count": 1024,
            "spectral_width": 20000.0,
            "reference_offset": 1220.7031,
            "origin_offset": 54237080.0,
        }
    ]
    experiment = cp.as_csdm(np.zeros(1024))
    sim.spin_systems.append(sys)
    sim.methods.append(
        BlochDecayCTSpectrum(
            channels=["2H"],
            rotor_frequency=5e3,
            spectral_dimensions=spectral_dims,
            experiment=experiment,
        )
    )
    sim.methods.append(
        BlochDecaySpectrum(
            channels=["2H"],
            rotor_frequency=10e3,
            spectral_dimensions=spectral_dims,
            experiment=experiment,
        )
    )
    return sim


def setup_signal_processor():
    op_list1 = [
        sp.IFFT(dim_index=0),
        sp.apodization.Exponential(FWHM="100 Hz"),
        sp.apodization.Gaussian(FWHM="200 Hz"),
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


def setup_second_signal_processor():
    op_list1 = [
        sp.IFFT(dim_index=0),
        sp.apodization.Exponential(FWHM="150 Hz"),
        sp.apodization.Gaussian(FWHM="250 Hz"),
        sp.FFT(dim_index=0),
        sp.Scale(factor=8),
        sp.baseline.ConstantOffset(offset=33.0),
        sp.Linear(amplitude=15.8, offset=11.1),
    ]
    op_list2 = [
        sp.Scale(factor=15),
        sp.baseline.ConstantOffset(offset=-33.0),
        sp.Linear(amplitude=1.1, offset=0.4),
    ]
    return [sp.SignalProcessor(operations=op) for op in [op_list1, op_list2]]


def setup_params():
    sim = setup_simulator()
    processor = setup_signal_processor()
    return sf.make_LMFIT_params(sim, processor, include="rotor_frequency")
