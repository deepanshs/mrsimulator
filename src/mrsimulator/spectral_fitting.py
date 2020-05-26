# -*- coding: utf-8 -*-
import numpy as np
from lmfit import Parameters
from mrsimulator import Simulator
from mrsimulator.post_simulation import Apodization
# from lmfit import minimize


__author__ = "Maxwell C Venetos"
__email__ = "maxvenetos@gmail.com"


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
    my_string = my_string.replace("isotopomers[", "ISO_")
    my_string = my_string.replace("].sites[", "_SITES_")
    my_string = my_string.replace(
        "].isotropic_chemical_shift", "_isotropic_chemical_shift"
    )
    my_string = my_string.replace("].shielding_symmetric.", "_shielding_symmetric_")
    my_string = my_string.replace("].quadrupolar.", "_quadrupolar_")
    my_string = my_string.replace("].abundance", "_abundance")
    my_string = my_string.replace("methods[", "METHODS_")  #
    my_string = my_string.replace("].post_simulation", "_POST_SIM_")
    my_string = my_string.replace(".scale", "scale")
    my_string = my_string.replace("['apodization'][", "APODIZATION_")
    my_string = my_string.replace("].args", "_args")

    return my_string


def _html_to_string(my_string):
    """
    Converts the HTML numbers to '[', ']', and '.' to allow for execution of the
    parameter name to update the simulator.

    Args:
        my_string: A string object

    Returns:
        String Object.

    """
    my_string = my_string.replace("ISO_", "isotopomers[")
    my_string = my_string.replace("_SITES_", "].sites[")
    my_string = my_string.replace(
        "_isotropic_chemical_shift", "].isotropic_chemical_shift"
    )
    my_string = my_string.replace("_shielding_symmetric_", "].shielding_symmetric.")
    my_string = my_string.replace("_quadrupolar_", "].quadrupolar.")
    my_string = my_string.replace("_abundance", "].abundance")

    my_string = my_string.replace("METHODS_", "methods[")  #
    my_string = my_string.replace("_POST_SIM_", "].post_simulation")
    my_string = my_string.replace("scale", "['scale']")
    my_string = my_string.replace("APODIZATION_", "['apodization'][")
    my_string = my_string.replace("_args", "].args")

    return my_string


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
    for i in range(len(sim.methods)):
        if sim.methods[i].post_simulation is not None:
            parent = f"methods[{i}].post_simulation"
            temp_list += [
                item
                for item in _traverse_dictionaries(
                    sim.methods[0].post_simulation, parent=parent
                )
                if "scale" in item
            ]
            if sim.methods[i].post_simulation["apodization"] is not None:
                for j in range(len(sim.methods[i].post_simulation["apodization"])):
                    temp_list.append(
                        _str_to_html(parent + f"['apodization'][{j}].args")
                    )

    length = len(sim.isotopomers)
    abundance = 0
    last_abund = f"{length - 1}_abundance"
    expression = "100"
    for i in range(length - 1):
        expression += f"-ISO_{i}_abundance"
    for i in range(length):
        abundance += eval("sim." + _html_to_string(f"isotopomers[{i}].abundance"))

    for items in temp_list:
        if "_eta" in items or "abundance" in items and last_abund not in items:
            if "_eta" in items:
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
            value = eval("sim." + _html_to_string(items))
            if type(value) == list:
                params.add(name=items, value=value[0])
            else:
                params.add(name=items, value=value)

    return params


def min_function(params, sim, apodization_function=None):
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
    # if not isinstance(data, cp.CSDM):
    #     raise ValueError(f"Expecting a `CSDM` object, found {type(data).__name__}.")
    if not isinstance(sim, Simulator):
        raise ValueError(f"Expecting a `Simulator` object, found {type(sim).__name__}.")


    # intensity_data = data.dependent_variables[0].components[0].real
    values = params.valuesdict()
    for items in values:
        if "args" not in items:
            nameString = "sim." + _html_to_string(items)
            executable = f"{nameString} = {values[items]}"
            exec(executable)
        else:
            nameString = "sim." + _html_to_string(items)
            executable = f"{nameString} = [{values[items]}]"
            exec(executable)

    sim.run()
    residual = np.asarray([])

    for i, method in enumerate(sim.methods):
        y_factored = method.apodize().real
        residual = np.append(
            residual,
            method.experiment.dependent_variables[0].components[0].real - y_factored,
        )

    return residual
