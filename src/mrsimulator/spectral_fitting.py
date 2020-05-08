# -*- coding: utf-8 -*-
import csdmpy as cp
from lmfit import minimize
from lmfit import Parameters
from mrsimulator import Simulator
from mrsimulator.apodization import Apodization


__author__ = "Maxwell C Venetos"
__email__ = "maxvenetos@gmail.com"


# def line_broadening(csdm_obj, sigma, broadType):
#     """
#     Applies appodization filter to the simulated spectrum using either Lorentzian
#     filtering or Gaussian filtering.

#     Args:
#     Returns:
#         Array of appodized intensities.
#     TimeDomain = ifft(ifftshift(amp))
#     if broadType == 0 and sigma != 0:
#         broadSignal = np.exp(-sigma * np.pi * np.abs(time))
#     # Gaussian broadening:
#     elif broadType == 1 and sigma != 0:
#         broadSignal = np.exp(-((time * sigma * np.pi) ** 2) * 2)
#     # No apodization
#     else:
#         broadSignal = 1

#     appodized = np.roll(TimeDomain * broadSignal, -int(x.count / 2))

#     return fftshift(fft(appodized)).real


def _str_to_html(my_string):
    """
    LMFIT Parameters class does not allow for names to include special characters.
    This function converts '[', ']', and '.' to their HTML numbers to comply with
    LMFIT.

    Args:
        my_string: A string object

    Returns:
        String object.

    """

    return my_string.replace("[", "91").replace("]", "93").replace(".", "46")


def _html_to_string(my_string):
    """
    Converts the HTML numbers to '[', ']', and '.' to allow for execution of the
    parameter name to update the simulator.

    Args:
        my_string: A string object

    Returns:
        String Object.

    """
    return my_string.replace("91", "[").replace("93", "]").replace("46", ".")


def _list_of_dictionaries(my_list):
    """
    Helper function for traverse_dictionaries function which will return a list of
    dictionaries.

    Args:
        my_list: A list object

    Returns:
        List Object.

    """
    return [item.dict() for item in my_list]


exclude = ["property_units", "isotope", "name", "description", "transition_pathways"]


def _traverse_dictionaries(dictionary, parent="isotopomers"):
    """
    Parses through the dictionary objects contained within the simulator object in
    order to return a list of all attributes that are populated.

    Args:
        dictionary: A dictionary or list object of the Isotopomer attributes from a
        simulation object
        parent: a string object used to create the addresses of the Isotopomer
        attributes.

    Returns:
        List Object.

    """
    name_list = []
    if isinstance(dictionary, dict):
        for key, vals in dictionary.items():
            if key not in exclude and vals is not None:
                if isinstance(vals, (dict, list)):
                    name_list += _traverse_dictionaries(
                        vals, _str_to_html(f"{parent}.{key}")
                    )
                else:
                    name_list += [_str_to_html(f"{parent}.{key}")]
    elif isinstance(dictionary, list):
        for i, items in enumerate(dictionary):
            name_list += _traverse_dictionaries(items, _str_to_html(f"{parent}[{i}]"))

    else:
        name_list += [_str_to_html(f"{parent}.{dictionary}")]

    return name_list


def make_fitting_parameters(sim, exclude_key=None):
    """
    Parses through the fitting parameter list to create LMFIT parameters used for
    fitting.

    Args:
        sim: a Simulator object.

    Returns:
        LMFIT Parameters object.

    """
    if not isinstance(sim, Simulator):
        raise ValueError(f"Expecting a `Simulator` object, found {type(sim).__name__}.")

    params = Parameters()
    temp_list = _traverse_dictionaries(_list_of_dictionaries(sim.isotopomers))
    length = len(sim.isotopomers)
    abundance = 0
    last_abund = f"{length - 1}9346abundance"
    expression = "100"
    for i in range(length - 1):
        expression += f"-isotopomers91{i}9346abundance"
    for i in range(length):
        abundance += eval("sim." + _html_to_string(f"isotopomers91{i}9346abundance"))

    for items in temp_list:
        if "46eta" in items or "abundance" in items and last_abund not in items:
            if "46eta" in items:
                params.add(
                    name=items,
                    value=eval("sim." + _html_to_string(items)),
                    min=0,
                    max=1,
                )
            if "abundance" in items:
                params.add(
                    name=items,
                    value=eval("sim." + _html_to_string(items)) / abundance * 100,
                    min=0,
                    max=100,
                )
        elif last_abund in items:
            params.add(
                name=items,
                value=eval("sim." + _html_to_string(items)),
                min=0,
                max=100,
                expr=expression,
            )
        else:
            params.add(name=items, value=eval("sim." + _html_to_string(items)))

    return params


function_mapping = {
    "Gaussian": Apodization.Gaussian,
    "Lorentzian": Apodization.Lorentzian,
}


def min_function(params, data, sim, apodization_function=None):
    """
    The simulation routine to establish how the parameters will update the simulation.

    Args:
        params: Parameters object containing parameters to vary during minimization.
        data: a CSDM object of the data to fit the simulation to.
        sim: Simulator object to be fit to data. Initialized with arbitrary fitting
        parameters.
        apodization_function: A string indicating the apodization function to use.
        Currently "Gaussian" and "Lorentzian" are supported.

    Returns:
        Array of the differences between the simulation and the experimental data.

    """
    if not isinstance(params, Parameters):
        raise ValueError(
            f"Expecting a `Parameters` object, found {type(params).__name__}."
        )
    if not isinstance(data, cp.CSDM):
        raise ValueError(f"Expecting a `CSDM` object, found {type(data).__name__}.")
    if not isinstance(sim, Simulator):
        raise ValueError(f"Expecting a `Simulator` object, found {type(sim).__name__}.")
    if not isinstance(apodization_function, str):
        raise ValueError(
            f"Expecting a `string` object, found {type(apodization_function).__name__}."
        )

    # intensity_data = data.dependent_variables[0].components[0].real
    values = params.valuesdict()
    for items in values:
        if items not in ["sigma", "factor"]:
            nameString = "sim." + _html_to_string(items)
            executable = f"{nameString} = {values[items]}"
            exec(executable)

    sim.run()
    y = sim.apodize(function_mapping[apodization_function], sigma=values["sigma"])

    y_factored = y * values["factor"]

    return (
        data.dependent_variables[0].components[0].real - y_factored
    )  # _factored#simulatedData.dependent_variables[0].components[0]


def spectral_fitting(experiment, sim, apodization_function, params):
    """
    Spectrum fitting routine to fit the mrsimulation to an experimental spectrum.
    Parameters may be provided or if not provided will be generated based on the
    simulation object passed through.

    Returns:
        CSDM object containing the experimental data and the simulated fit.

    """
    intensity_data = experiment.dependent_variables[0].components[0].real
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
        min_function,
        params,
        fcn_args=(
            experiment.dependent_variables[0].components[0].real,
            sim,
            apodization_function,
        ),
    )
    result = minner.minimize()

    # report_fit(result)
    values = result.params

    sim.run()
    # sim_data = sim.as_csdm_object()

    if apodization_function is not None:
        y = sim.apodize(function_mapping[apodization_function], sigma=values["sigma"])
    else:
        y = sim.methods[0].simulation.to_list()[1]
        # y = sim.dependent_variables[0].components[0]

    if "factor" in values:
        y_factored = y * values["factor"]
    else:
        y_factored = y

    return (
        intensity_data - y_factored
    )  # _factored#simulatedData.dependent_variables[0].components[0]


# def spectral_fitting(experiment, sim, apodization_function, params):
#     """
#     Spectrum fitting routine to fit the mrsimulation to an experimental spectrum
# . Parameters may be provided or if not provided will be generated based on the
#  simulation object passed through.

#     Returns:
#         CSDM object containing the experimental data and the simulated fit.

#     """
#     if len(params) == 0:
#         params = make_fitting_parameters(sim)
#         params.add(
#             name="sigma", value=experiment.dimensions[0].
# increment.to("Hz").value, min=0
#         )
#         params.add(
#             name="factor",
#             value=experiment.dependent_variables[0].components[0].max().real,
#             min=0,
#         )
#     if "sigma" not in params:
#         params.add(
#             name="sigma", value=experiment.dimensions[0].increment.to("Hz").value,
#  min=0
#         )
#     if "factor" not in params:
#         params.add(
#             name="factor",
#             value=experiment.dependent_variables[0].components[0].max().real,
#             min=0,
#         )

#     minner = Minimizer(
#         min_function,
#         params,
#         fcn_args=(
#             experiment.dependent_variables[0].components[0].real,
#             sim,
#             apodization_function,
#         ),
#     )
#     result = minner.minimize()

#     report_fit(result)

#     sim.run(method=one_d_spectrum)
#     sim_data = sim.as_csdm_object()

#     x_simulated = sim_data.dimensions[0]
#     y_simulated = sim_data.dependent_variables[0].components[0]
#     sim_data.dependent_variables[0].components[0] = line_broadening(
#         x_simulated, y_simulated, result.params["sigma"], 0
#     ).real
#     sim_data.dependent_variables[0].components[0] = (
#         sim_data.dependent_variables[0].components[0] * result.params["factor"]
#     )

#     new_csdm = cp.new()

#     new_csdm.add_dimension(experiment.dimensions[0])
#     new_csdm.add_dependent_variable(experiment.dependent_variables[0])

#     new_csdm.add_dependent_variable(sim_data.dependent_variables[0])

#     return new_csdm
