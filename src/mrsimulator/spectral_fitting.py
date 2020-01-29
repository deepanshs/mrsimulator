# -*- coding: utf-8 -*-
from pprint import pprint

import csdmpy as cp
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from lmfit import Minimizer
from lmfit import Parameters
from lmfit import report_fit
from mrsimulator import Dimension
from mrsimulator import Isotopomer
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SymmetricTensor as st
from mrsimulator.apodization import Apodization
from mrsimulator.methods import one_d_spectrum
from numpy.fft import fft
from numpy.fft import fftshift
from numpy.fft import ifft
from numpy.fft import ifftshift


__author__ = "Maxwell C Venetos"
__email__ = "maxvenetos@gmail.com"


def line_broadening(x, amp, sigma, broadType):
    """
    Applies appodization filter to the simulated spectrum using either Lorentzian filtering or Gaussian filtering.

    Returns:
        Array of appodized intensities.

    """
    freq = x.coordinates.to("Hz")

    TimeDomain = ifft(ifftshift(amp))
    TimeDomain = np.roll(TimeDomain, int(x.count / 2))
    t = np.arange(x.count) - int(x.count / 2)

    time = t * 1 / (len(freq) * x.increment.to("Hz").value)

    # Lorentzian broadening:
    if broadType == 0 and sigma != 0:
        broadSignal = np.exp(-sigma * np.pi * np.abs(time))
    # Gaussian broadening:
    elif broadType == 1 and sigma != 0:
        broadSignal = np.exp(-((time * sigma * np.pi) ** 2) * 2)
    # No apodization
    else:
        broadSignal = 1

    appodized = np.roll(TimeDomain * broadSignal, -int(x.count / 2))

    return fftshift(fft(appodized))


def str_to_html(my_string):
    """
    LMFIT Parameters class does not allow for names to include special characters. This function converts '[', ']', and '.' to their HTML numbers to comply with LMFIT.

    Returns:
        String object.

    """

    return my_string.replace("[", "91").replace("]", "93").replace(".", "46")


def html_to_string(my_string):
    """
    Converts the HTML numbers to '[', ']', and '.' to allow for execution of the parameter name to update the simulator.

    Returns:
        String Object.

    """
    return my_string.replace("91", "[").replace("93", "]").replace("46", ".")


def list_of_dictionaries(my_list):
    """
    Helper function for traverse_dictionaries function which will return a list of dictionaries.

    Returns:
        List Object.

    """
    return [item.dict() for item in my_list]


exclude = ["property_units", "isotope", "name", "description"]


def traverse_dictionaries(dictionary, parent="isotopomers"):
    """
    Parses through the dictionary objects contained within the simulator object in order to return a list of all attributes that are populated.

    Returns:
        List Object.

    """
    name_list = []
    if isinstance(dictionary, dict):
        for key, vals in dictionary.items():
            if key not in exclude and vals is not None:
                if isinstance(vals, (dict, list)):
                    name_list += traverse_dictionaries(
                        vals, str_to_html(f"{parent}.{key}")
                    )
                else:
                    name_list += [str_to_html(f"{parent}.{key}")]
    elif isinstance(dictionary, list):
        for i, items in enumerate(dictionary):
            name_list += traverse_dictionaries(items, str_to_html(f"{parent}[{i}]"))

    else:
        name_list += [str_to_html(f"{parent}.{dictionary}")]

    return name_list


def make_fitting_parameters(sim):
    """
    Parses through the fitting parameter list to create LMFIT parameters used for fitting.

    Returns:
        LMFIT Parameters object.

    """
    params = Parameters()
    temp_list = traverse_dictionaries(list_of_dictionaries(sim.isotopomers))
    length = len(sim.isotopomers)
    abundance = 0
    last_abund = f"{length - 1}9346abundance"
    expression = "100"
    for i in range(length - 1):
        expression += f"-isotopomers91{i}9346abundance"
    for i in range(length):
        abundance += eval("sim." + html_to_string(f"isotopomers91{i}9346abundance"))

    for items in temp_list:
        if "46eta" in items or "abundance" in items and last_abund not in items:
            if "46eta" in items:
                params.add(
                    name=items, value=eval("sim." + html_to_string(items)), min=0, max=1
                )
            if "abundance" in items:
                params.add(
                    name=items,
                    value=eval("sim." + html_to_string(items)) / abundance * 100,
                    min=0,
                    max=100,
                )
        elif last_abund in items:
            params.add(
                name=items,
                value=eval("sim." + html_to_string(items)),
                min=0,
                max=100,
                expr=expression,
            )
        else:
            params.add(name=items, value=eval("sim." + html_to_string(items)))

    return params


function_mapping = {
    "Gaussian": Apodization.Gaussian,
    "Lorentzian": Apodization.Lorentzian,
}


def fcn2min(params, data, sim, apodization_function):
    """
    The simulation routine to establish how the parameters will update the simulation.

    Returns:
        Array of the differences between the simulation and the experimental data.

    """
    values = params.valuesdict()
    for items in values:
        if items not in ["sigma", "factor"]:
            nameString = "sim." + html_to_string(items)
            executable = f"{nameString} = {values[items]}"
            exec(executable)

    sim.run(method=one_d_spectrum)
    y = sim.apodize(function_mapping[apodization_function], sigma=values["sigma"])

    y_factored = y * values["factor"]

    return data - y_factored  # _factored#simulatedData.dependent_variables[0].components[0]


def spectral_fitting(experiment, sim, apodization_function, params):
    """
    Spectrum fitting routine to fit the mrsimulation to an experimental spectrum. Parameters may be provided or if not provided will be generated based on the simulation object passed through.

    Returns:
        CSDM object containing the experimental data and the simulated fit.

    """
    if len(params) == 0:
        params = make_fitting_parameters(sim)
        params.add(
            name="sigma", value=experiment.dimensions[0].increment.to("Hz").value, min=0
        )
        params.add(
            name="factor",
            value=experiment.dependent_variables[0].components[0].max().real,
            min=0,
        )
    if "sigma" not in params:
        params.add(
            name="sigma", value=experiment.dimensions[0].increment.to("Hz").value, min=0
        )
    if "factor" not in params:
        params.add(
            name="factor",
            value=experiment.dependent_variables[0].components[0].max().real,
            min=0,
        )

    minner = Minimizer(
        fcn2min,
        params,
        fcn_args=(
            experiment.dependent_variables[0].components[0].real,
            sim,
            apodization_function,
        ),
    )
    result = minner.minimize()

    report_fit(result)

    sim.run(method=one_d_spectrum)
    sim_data = sim.as_csdm_object()

    x_simulated = sim_data.dimensions[0]
    y_simulated = sim_data.dependent_variables[0].components[0]
    sim_data.dependent_variables[0].components[0] = line_broadening(
        x_simulated, y_simulated, result.params["sigma"], 0
    ).real
    sim_data.dependent_variables[0].components[0] = (
        sim_data.dependent_variables[0].components[0] * result.params["factor"]
    )

    new_csdm = cp.new()

    new_csdm.add_dimension(experiment.dimensions[0])
    new_csdm.add_dependent_variable(experiment.dependent_variables[0])

    new_csdm.add_dependent_variable(sim_data.dependent_variables[0])

    return new_csdm
