# -*- coding: utf-8 -*-
import mrsimulator.signal_processing as sp
import mrsimulator.signal_processing.apodization as apo
from mrsimulator import Simulator

try:
    from lmfit import Parameters

    FOUND_LMFIT = True
except ImportError:
    FOUND_LMFIT = False


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
    my_string = my_string.replace("spin_systems[", "ISO_")
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
    my_string = my_string.replace(".apodization[", "APODIZATION_")
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
    my_string = my_string.replace("ISO_", "spin_systems[")
    my_string = my_string.replace("_SITES_", "].sites[")
    my_string = my_string.replace(
        "_isotropic_chemical_shift", "].isotropic_chemical_shift"
    )
    my_string = my_string.replace("_shielding_symmetric_", "].shielding_symmetric.")
    my_string = my_string.replace("_quadrupolar_", "].quadrupolar.")
    my_string = my_string.replace("_abundance", "].abundance")

    my_string = my_string.replace("METHODS_", "methods[")  #
    my_string = my_string.replace("_POST_SIM_", "].post_simulation")
    my_string = my_string.replace("scale", ".scale")
    my_string = my_string.replace("APODIZATION_", ".apodization[")
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


exclude = [
    "property_units",
    "isotope",
    "name",
    "label",
    "description",
    "transition_pathways",
]


def _traverse_dictionaries(dictionary, parent="spin_systems"):
    """
    Parses through the dictionary objects contained within the simulator object in
    order to return a list of all attributes that are populated.

    Args:
        dictionary: A dictionary or list object of the SpinSystem attributes from a
        simulation object
        parent: a string object used to create the addresses of the SpinSystem
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

    return name_list


def _post_sim_LMFIT_params(post_sim):
    """
    Creates an LMFIT Parameters object for SignalProcessor operations
    involved in spectrum fitting

    post_sim: SignalProcessor object

    returns: Parameters object
    """
    temp_dict = {}
    # for item in post_sim.operations:
    #     prepend = f"DEP_VAR_{item.dependent_variable}_"
    for i, operation in enumerate(post_sim.operations):
        if isinstance(operation, apo.Gaussian):
            identifier = f"operation_{i}_Gaussian"
            arg = operation.sigma
            temp_dict[f"{identifier}"] = arg
        elif isinstance(operation, apo.Exponential):
            identifier = f"operation_{i}_Exponential"
            arg = operation.Lambda
            temp_dict[f"{identifier}"] = arg
        elif isinstance(operation, sp.Scale):
            identifier = f"operation_{i}_Scale"
            arg = operation.factor
            temp_dict[f"{identifier}"] = arg

    params = Parameters()
    for key, val in temp_dict.items():
        params.add(name=key, value=val)

    return params


def _update_post_sim_from_LMFIT_params(params, post_sim):
    """
    Updates SignalProcessor operation arguments from an
    LMFIT Parameters object

    params: LMFIT Parameters object
    post_sim: SignalProcessor object

    """
    temp_dict = {}
    arg_dict = {"Gaussian": "sigma", "Exponential": "Lambda", "Scale": "factor"}
    for param in params:
        # iterating through the parameter list looking for only DEP_VAR (ie post_sim params)
        if "operation_" in param:
            # splitting parameter name to obtain
            # Dependent variable index (var)
            # index of operation in the operation list (opIndex)
            # arg value for the operation (val)
            split_name = param.split("_")
            # var = split_name[split_name.index("VAR") + 1]
            opIndex = split_name[split_name.index("operation") + 1]
            val = params[param].value
            # creating a dictionary of operations and arguments for each dependent variable
            # if f"DepVar_{var}" not in temp_dict.keys():
            #     temp_dict[f"DepVar_{var}"] = {}
            temp_dict[f"{opIndex}_{split_name[-1]}"] = val

    # iterate through list of operation lists
    # for item in post_sim.operations:
    # iterating through dictionary with corresponding dependent variable index
    for operation, val in temp_dict.items():
        # creating assignment strings to create the correct address for updating each operation
        split = operation.split("_")
        # dep_var_operation_list = f"post_sim.operations[{item.dependent_variable}]"
        operation_val_update = f"post_sim.operations[{split[0]}].{arg_dict[split[-1]]}"
        assignment = f"={val}"
        exec(operation_val_update + assignment)


def make_LMFIT_parameters(sim, post_sim=None, exclude_key=None):
    """
    Parses through the fitting parameter list to create LMFIT parameters used for
    fitting.

    Args:
        sim: a Simulator object.
        post_sim: a SignalProcessor object

    Returns:
        LMFIT Parameters object.

    """
    if not FOUND_LMFIT:
        error = (
            f"The helper function {__name__} requires 'lmfit' module to create lmfit "
            "paramters. Please install the lmfit module using\n'pip install lmfit'.",
        )
        raise ImportError(error)

    if not isinstance(sim, Simulator):
        raise ValueError(f"Expecting a `Simulator` object, found {type(sim).__name__}.")
    if not isinstance(post_sim, sp.SignalProcessor) or post_sim is None:
        raise ValueError(
            f"Expecting a `SignalProcessor` object, found {type(post_sim).__name__}."
        )

    if isinstance(sim, Simulator):
        params = Parameters()
        temp_list = _traverse_dictionaries(_list_of_dictionaries(sim.spin_systems))

        length = len(sim.spin_systems)
        abundance = 0
        last_abund = f"{length - 1}_abundance"
        expression = "100"
        for i in range(length - 1):
            expression += f"-ISO_{i}_abundance"
        for i in range(length):
            abundance += eval("sim." + _html_to_string(f"spin_systems[{i}].abundance"))

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
    else:
        params = Parameters()
        temp_list = _traverse_dictionaries(_list_of_dictionaries(sim))

    if isinstance(post_sim, sp.SignalProcessor):
        temp_params = _post_sim_LMFIT_params(post_sim)
        for item in temp_params:
            params.add(name=item, value=temp_params[item].value)
        # params.add_many(temp_params)

    return params


def LMFIT_min_function(params, sim, post_sim=None):
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
    if not isinstance(post_sim, sp.SignalProcessor) or post_sim is None:
        raise ValueError(
            f"Expecting a `SignalProcessor` object, found {type(post_sim).__name__}."
        )

    if not isinstance(sim, Simulator):
        raise ValueError(f"Expecting a `Simulator` object, found {type(sim).__name__}.")

    values = params.valuesdict()
    for items in values:
        if "operation_" not in items:
            nameString = "sim." + _html_to_string(items)
            executable = f"{nameString} = {values[items]}"
            exec(executable)
        elif "operation_" in items and post_sim is not None:
            _update_post_sim_from_LMFIT_params(params, post_sim)

    sim.run()
    post_sim.data = sim.methods[0].simulation
    processed_data = post_sim.apply_operations()
    # residual = np.asarray([])

    if sim.config.decompose_spectrum == "spin_system":
        datum = 0
        for decomposed_datum in processed_data.dependent_variables:
            datum += decomposed_datum.components[0]
            # datum = [sum(i) for i in zip(datum, decomposed_datum)]
    else:
        datum = processed_data.dependent_variables[0].components[0]

    return (
        sim.methods[0].experiment.dependent_variables[0].components[0].real - datum.real
    )

    # MULTIPLE EXPERIMENTS
    # for i, method in enumerate(sim.methods):
    #     y_factored = method.apodize().real
    #     residual = np.append(
    #         residual,
    #         method.experiment.dependent_variables[0].components[0].real - y_factored,
    #     )

    # return residual
