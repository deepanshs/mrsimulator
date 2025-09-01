# -*- coding: utf-8 -*-
import copy
import re

import emcee
import numpy as np
from mrsimulator.utils.spectral_fitting import get_correct_data_order

__author__ = ["He Sun", "Deepansh Srivastava"]
__email__ = ["wushanyun64@gmail.com", "srivastava.89@osu.edu"]


def name_abbrev(params):
    """
    The limfit parameters obj created by function make_LMFIT_params contains detailed
    but lengthy names which could be hard for ploting.
    This small function simplifys the long default param name for substream ploting.

    Example (default name --> abbreviated name)
    -----------------------------
    sys_0_site_0_isotropic_chemical_shift --> delta_0_0
    sys_1_site_0_shielding_symmetric_eta-->etaCS_1_0
    mth_0_rotor_frequency-->spin_freq_0
    SP_0_operation_2_Exponential_FWHM -->expo_fwhm_0_2

    Paramsters
    -----------------------------
    params: LMFIT Parameters.
        The lmfit parameters obj with default mrsimulator param name.

    Return
    -----------------------------
    name_abbrev_list: list
        A list of abbreviated parameters names.
    """
    abbreviation_pairs = {
        r"sys_[0-9]+_site_[0-9]+_isotropic_chemical_shift": "delta_iso",
        r"sys_[0-9]+_site_[0-9]+_shielding_symmetric_zeta": "zeta",
        r"sys_[0-9]+_site_[0-9]+_shielding_symmetric_eta": "etaCS",
        r"sys_[0-9]+_site_[0-9]+_quadrupolar_Cq": "Cq",
        r"sys_[0-9]+_site_[0-9]+_quadrupolar_eta": "etaQ",
        r"mth_[0-9]+_rotor_frequency": "spin_freq",
        r"sys_[0-9]+_abundance": "abundance",
        r"SP_[0-9]+_operation_[0-9]+_Exponential_FWHM": "expo_fwhm",
        r"SP_[0-9]+_operation_[0-9]+_Gaussian_FWHM": "gauss_fwhm",
        r"SP_[0-9]+_operation_[0-9]+_Scale_factor": "scale",
        r"SP_[0-9]+_operation_[0-9]+_ConstantOffset_offset": "offset",
        r"SP_[0-9]+_operation_[0-9]+_Linear_amplitude": "linear",
    }

    name_abbrev_list = []
    for name in params.keys():
        for pattern, abbrev in abbreviation_pairs.items():
            pattern = re.compile(pattern)
            if pattern.fullmatch(name):
                numbers = re.findall(r"(_[0-9]+)_", name)
                name_abbrev = abbrev
                for num in numbers:
                    name_abbrev += num
                name_abbrev_list.append(name_abbrev)
    return name_abbrev_list


class mrsim_emcee:
    """
    A utility class to sample the posterior distribution function (PDF) of NMR
    paramters with Markov Chain Monte Carlo (MCMC) method.

    This is a class depends on the emcee package.
    https://emcee.readthedocs.io/en/stable/#

    For more information on Markov Chain Monte Carlo, check:
    https://knac.middlebury.edu/symp_sessions/2015/flaherty_MCMC_2015.pdf

    Parameters
    --------------------------
    params: LMFIT parameters
        Initial guess of the parameters to start MCMC walkers. User can use
        mrsimulator.utils.spectral_fitting.make_LMFIT_params to generate LMFIT
        parameters for mrsimulator. Please note mrsim_emcee is not a fitting
        package, the aim here to sample the posterior distribution of the NMR
        parameters AFTER fitting. We recommend user to first get a proper fitting
        with spectral_fitting model then start this process with the fitted params.
    simulator: MRsimulator simulator
        The simulator object constructed with NMR paramters.
    processor: MRsimulator processor
        The processor object storing NMR signal processing parameters.
    sigma: float
        The noise level of the spectrum baseline.
    """

    def __init__(self, params, simulator, processor, sigma=None):
        params = copy.deepcopy(params)
        self.params = self._vary_only(params)
        self.nvarys = len(self.params)
        self.simulator = simulator
        self.processor = processor
        self.sigma = sigma

    @staticmethod
    def _vary_only(params):
        """
        Filter the parameter list, leave only variables.

        Note: abundance should not a a variable, it will be removed from the params
        even if included by user with Vary = True.

        Parameters
        ------------------------------
        params: Parameters
            LMFIT parameters obj that store the variables.
        """
        for k, v in list(params.items()):
            if not v.vary or "abundance" in k:
                params.pop(k)
        return params

    def mcmc(self, steps=1000, nwalkers=100, burn=100, thin=10, progress=True):
        """
        Bayesian sampling of the posterior distribution.

        This method uses the emcee package 'https://emcee.readthedocs.io/en/stable/'
        to sample the posterior distribution of NMR parameters with Marcov Chain
        Monte Carlo (MCMC) method.

        Parameters
        -----------------------
        steps: int
            The total number of samples to draw from the posterior distribution
            for each walkers.
        nwalkers: int
            The numbers of walkers in the emnsemble.
        burn:
            Number of steps to discard at the begining of the sampling process.
            The first few steps when the chain is making its way towards the
            center of the PDF, are called ‘burn in’ and are removed since they
            do not represent the final probability density function.
        thin: int
            Only accept 1 in every `thin` samples.
        progress: bool
            If True, show a progress bar while running.

        Return
        ---------------------------
        Result: dictionary
            Result have 4 keys.
            'chain' contains the original chain of the random walker for each paramters
            and walkers, dimention [(steps - burn) // thin, nwalkers, nvarys].
            'flat_chain' contians the faltened chain the chain for all the walkers
            on each parameters.
            'log_prob' contains output _log_probability function.
            'variables' contains all the variables' (none variables not included)
            values, std, relative error, initial values etc.
        """
        # randstate
        seed = 100
        rand = np.random.RandomState(seed)
        # initial guess
        var = np.zeros(self.nvarys)
        for i, v in enumerate(self.params):
            param = self.params[v]
            var[i] = param.value
        p0 = 1 + rand.randn(nwalkers, self.nvarys) * 1.0e-4
        p0 *= var

        # run mcmc sampler
        sampler = emcee.EnsembleSampler(
            nwalkers,
            self.nvarys,
            self._log_probability,
            args=(self.params, self.simulator, self.processor, self.sigma),
        )

        sampler.run_mcmc(p0, steps, progress=progress)

        result = {}
        result["chain"] = sampler.get_chain(thin=thin, discard=burn, flat=False)
        result["flat_chain"] = sampler.get_chain(thin=thin, discard=burn, flat=True)
        result["log_prob"] = sampler.get_log_prob(thin=thin, discard=burn)
        result["accept_frac"] = sampler.acceptance_fraction
        quantiles = np.percentile(result["flat_chain"], [15.87, 50, 84.13], axis=0)
        result["variables"] = self.params

        for i, k in enumerate(result["variables"]):
            std_l, median, std_u = quantiles[:, i]
            result["variables"][k].value = median
            result["variables"][k].stderr = 0.5 * (std_u - std_l)

        return result

    @staticmethod
    def _log_prior(theta, bounds):
        """
        Calculate the log prior for the parameters given.
        If the values in params are out side its bound, return -inf, else return 0.

        Parameters
        -----------------------
        theta: sequence
            sequence of the variables.
        bounds: array
            Maximum and minimum of the variables.
        """
        if np.any(theta > bounds[0, :]) or np.any(theta < bounds[1, :]):
            return -np.inf
        return 0.0

    def _log_probability(self, theta, params, simulator, processor, sigma):
        """
        Calculate the log posterior probability function of the nmr parameters.
        """
        # calculate log prior.
        max_ = [params[k].max for k in params]
        min_ = [params[k].min for k in params]
        bounds = np.array([max_, min_])
        lp = self._log_prior(theta, bounds)
        if not np.isfinite(lp):
            return -np.inf

        for k, v in zip(params.keys(), theta):
            params[k].value = v
        fxn_output = self.minimization_function(params, simulator, processor, sigma)
        lnprob = np.asarray(fxn_output).ravel()
        lnprob = -0.5 * (lnprob * lnprob).sum()
        return lnprob

    @staticmethod
    def _update_sim_params(simulator, params):
        """
        update the simulation object with the NMR parameters provided.
        """
        simulator = copy.deepcopy(simulator)
        values = params.valuesdict()

        NMR_params = {
            "shielding_symmetric": ["zeta", "eta"],
            "quadrupolar": ["Cq", "eta", "beta"],
        }

        # Update simulation parameters iso, eta, and zeta for the site object
        num_spin = len(simulator.spin_systems)
        for i in range(num_spin):
            num_site = len(simulator.spin_systems[i].sites)
            for j in range(num_site):
                site = simulator.spin_systems[i].sites[j]
                # Update NMR params
                for k, v in NMR_params.items():
                    if getattr(site, k):
                        for p in v:
                            name = f"sys_{i}_site_{j}_{k}_{p}"
                            if getattr(getattr(site, k), p) and name in values.keys():
                                setattr(getattr(site, k), p, values[name])
                # Update isotropic chemical shift
                name = f"sys_{i}_site_{j}_isotropic_chemical_shift"
                if site.isotropic_chemical_shift and name in values.keys():
                    site.isotropic_chemical_shift = values[name]
        return simulator

    @staticmethod
    def _update_methods(simulator, params):
        """
        update the simulation object with the method parameters provided.
        """
        simulator = copy.deepcopy(simulator)
        values = params.valuesdict()
        # Update the spinning freq here.
        for i, method in enumerate(simulator.methods):
            for sd in method.spectral_dimensions:
                for e in sd.events:
                    if f"mth_{i}_rotor_frequency" in values:
                        e.rotor_frequency = values[f"mth_{i}_rotor_frequency"]
        return simulator

    @staticmethod
    def _update_signal_processors(processors, params):
        """
        update the MP signal processor object with the parameters provided.
        """
        values = params.valuesdict()

        # Define the dict of possible signal processors.
        processors_dict = {
            "Gaussian": "FWHM",
            "Exponential": "FWHM",
            "Scale": "factor",
            "ConstantOffset": "offset",
            "Linear": ["amplitude", "offset"],
        }

        # update the SignalProcessor parameter and apply line broadening.
        if isinstance(processors, list):
            processors = processors
        else:
            processors = [processors]

        for i, proc in enumerate(processors):
            for j, oper in enumerate(proc.operations):
                oper_name = oper.__class__.__name__
                if oper_name in processors_dict:
                    if oper_name == "Linear":
                        for attr in processors_dict[oper_name]:
                            setattr(
                                oper,
                                attr,
                                values[f"SP_{i}_operation_{j}_{oper_name}_{attr}"],
                            )
                    else:
                        attr = processors_dict[oper_name]
                        setattr(
                            oper,
                            attr,
                            values[f"SP_{i}_operation_{j}_{oper_name}_{attr}"],
                        )
        return processors

    def minimization_function(self, params, simulator, processors, sigma):
        """
        Definition of the minimization function used to fit the experimental spectrum
        """

        # Update simulation parameters
        simulator = self._update_sim_params(simulator, params)
        simulator = self._update_methods(simulator, params)

        # run the simulation
        simulator.run()

        # update processors
        processors = self._update_signal_processors(processors, params)

        # apply signal processing
        processed_data = [
            proc.apply_operations(method.simulation)
            for proc, method in zip(processors, simulator.methods)
        ]
        # set sigma as 1 if not defined
        sigma = [1.0 for _ in simulator.methods] if sigma is None else sigma
        sigma = sigma if isinstance(sigma, list) else [sigma]

        # calculate the residual.
        diff = np.asarray([])
        for processed_datum, mth, sigma_ in zip(
            processed_data, simulator.methods, sigma
        ):
            datum = 0
            for decomposed_datum in processed_datum.y:
                datum += decomposed_datum.components[0].real

            exp_data = get_correct_data_order(mth)
            diff = np.append(diff, (exp_data - datum) / sigma_)
        return diff
