import os
import csdmpy as cp
from mrsimulator.utils import get_spectral_dimensions
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import emcee
from mrsimulator import Simulator, SpinSystem, Site
from mrsimulator.methods import BlochDecaySpectrum
from mrsimulator import signal_processing as sp
from lmfit import Minimizer, Parameters, fit_report, minimize
import json
from mrsimulator.utils.spectral_fitting import get_correct_data_order
from multiprocessing import Pool
import copy

class mrsim_emcee:
    '''
    A utility class to sample the posterior distribution function (PDF) of NMR paramters with Markov Chain Monte Carlo (MCMC) method.
    This is a class depends on the emcee (https://emcee.readthedocs.io/en/stable/#) package.

    Parameters
    --------------------------
    params: LMFIT parameters
        Initial guess of the parameters to start MCMC walkers. User can use mrsimulator.utils.spectral_fitting.make_LMFIT_params to generate LMFIT parameters for mrsimulator.
    simulator: MRsimulator simulator
        The simulator object constructed with NMR paramters.
    processor: MRsimulator processor
        The processor object storing NMR signal processing parameters.
    sigma: float
        The noise level of the spectrum baseline.

    '''
    def __init__(self,params, simulator, processor, sigma = None):
        self.params = self._vary_only(params)
        self.nvarys = len(self.params)
        self.simulator = simulator
        self.processor = processor
        self.sigma = sigma

    @staticmethod
    def _vary_only(params):
        '''
        Filter the parameter list, leave only variables.

        Parameters
        ------------------------------
        params: Parameters
            LMFIT parameters obj that store the variables.
        '''
        for k, v in list(params.items()):
            if not v.vary:
                params.pop(k)
        return(params)

    def mcmc(self,steps = 1000, nwalkers = 100, burn = 100, thin = 10, progress = True):
        '''
        Bayesian sampling of the posterior distribution.

        This method uses the emcee package 'https://emcee.readthedocs.io/en/stable/' to sample the posterior distribution of NMR parameters with Marcov Chain Monte Carlo (MCMC) method.

        Parameters
        -----------------------
        steps: int
            The total number of samples to draw from the posterior distribution for each walkers.
        nwalkers: int
            The numbers of walkers in the emnsemble.
        burn:
            Number of steps to discard at the begining of the sampling process. The first few steps when the chain is making its way towards the center of the PDF,
            are called ‘burn in’ and are removed since they do not represent the final PDF.
        thin: int
            Only accept 1 in every `thin` samples.
        progress: bool
            If True, show a progress bar while running.

        '''
        # randstate
        seed = 100
        rand = np.random.RandomState(seed)
        # initial guess
        var = np.zeros(self.nvarys)
        for i, v in enumerate(self.params):
            param = self.params[v]
            var[i] = param.value
        p0 = 1 + rand.randn(nwalkers, self.nvarys) * 1.e-4
        p0*=var

        #run mcmc sampler
        with Pool() as pool:
            self.sampler = emcee.EnsembleSampler(nwalkers,self.nvarys,self._log_probability,
                                                 args = (self.params, self.simulator, self.processor, self.sigma),
                                                 pool = pool)

            self.sampler.run_mcmc(p0, steps, progress=progress)

        result = []
        result['flat_chain'] = self.sampler.get_chain(thin=thin, discard=burn, flat = True)
        result['log_prob'] = self.sampler.get_log_prob(thin=thin, discard = burn)
        quantiles = np.percentile(result['flat_chain'],[15.87,50,84.13],axis=0)
        result['params'] = copy.deepcopy(self.params)

        for i,k in enumerate(result['params']):
            std_l, median, std_u = quantiles[:,i]
            result['params'][k].value = median
            result['params'][k].stderr = 0.5*(std_u-std_l)

        return(result)

    @staticmethod
    def _log_prior(theta,bounds):
        '''
        Calculate the log prior for the parameters given.
        If the values in params are out side its bound, return -inf, else return 0.
        '''
        if np.any(theta>bounds[0,:]) or np.any(theta<bounds[1,:]):
                return -np.inf
        return 0.0

    def _log_probability(self,theta,params,simulator,processor,sigma):
        '''
        Calculate the log posterior probability function of the nmr parameters.

        '''
        #calculate log prior.
        max_ = [params[k].max for k in params]
        min_ = [params[k].min for k in params]
        bounds = np.array([max_,min_])
        lp = self._log_prior(theta,bounds)
        if not np.isfinite(lp):
            return -np.inf

        for k,v in zip(params.keys(),theta):
            params[k].value = v

        fxn_output = self.minimization_function(params, sim, processor, sigma)
        lnprob = np.asarray(fxn_output).ravel()
        lnprob = -0.5*(lnprob*lnprob).sum()
        return lnprob

    @staticmethod
    def minimization_function(params, simulator, processors, sigma=None):
        '''
        Definition of the minimization function used to fit the experimental spectrum
        '''
        values = params.valuesdict()

        # Update simulation parameters iso, eta, and zeta for the site object
        num_spin = len(simulator.spin_systems)
        for i in range(num_spin):
            num_site = len(simulator.spin_systems[i].sites)
            for j in range(num_site):
                site = simulator.spin_systems[i].sites[j]
                site.isotropic_chemical_shift = values[f"sys_{i}_site_{j}_isotropic_chemical_shift"]
                site.shielding_symmetric.eta = values[f"sys_{i}_site_{j}_shielding_symmetric_eta"]
                site.shielding_symmetric.zeta = values[f"sys_{i}_site_{j}_shielding_symmetric_zeta"]

        # Update the spinning freq here.
        for i, method in enumerate(simulator.methods):
            for sd in method.spectral_dimensions:
                for e in sd.events:
                    if f'mth_{i}_rotor_frequency' in values:
                        e.rotor_frequency = values[f'mth_{i}_rotor_frequency']

        # run the simulation
        simulator.run()

        # Define the dict of possible signal processors.
        processors_dict = {
            "Gaussian": "FWHM",
            "Exponential": "FWHM",
            "Scale": "factor",
        }

        # update the SignalProcessor parameter and apply line broadening.
        if isinstance(processors,list):
            processors = processors
        else:
            processors = [processors]

        for i, proc in enumerate(processors):
            for j, oper in enumerate(proc.operations):
                if oper.__class__.__name__ in processors_dict:
                    attr = processors_dict[oper.__class__.__name__]
                    setattr(oper,attr,values[f'SP_{i}_operation_{j}_{oper.__class__.__name__}_{attr}'])

        # apply signal processing
        processed_data = [proc.apply_operations(method.simulation) for proc, method in zip(processors, simulator.methods)]

        sigma = [1.0 for _ in sim.methods] if sigma is None else sigma
        sigma = sigma if isinstance(sigma, list) else [sigma]

        diff = np.asarray([])
        for processed_datum, mth, sigma_ in zip(processed_data, simulator.methods, sigma):
            datum = 0
            for decomposed_datum in processed_datum.y:
                datum += decomposed_datum.components[0].real

            exp_data = get_correct_data_order(mth)
            diff = np.append(diff, (exp_data - datum) / sigma_)
        return diff
