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

START = "sys_"
ENCODING_PAIRS = [
    ["spin_systems[", START],
    ["].sites[", "_site_"],
    ["].isotropic_chemical_shift", "_isotropic_chemical_shift"],
    ["].shielding_symmetric.", "_shielding_symmetric_"],
    ["].quadrupolar.", "_quadrupolar_"],
    ["].abundance", "_abundance"],
    ["methods[", "METHODS_"],  # why does methods need to be parameterized?
    ["].post_simulation", "_POST_SIM_"],
    [".scale", "scale"],
    [".apodization[", "APODIZATION_"],
    ["].args", "_args"],
]

EXCLUDE = [
    "property_units",
    "isotope",
    "name",
    "label",
    "description",
    "transition_pathways",
]


def _str_encode(my_string):
    """
    LMFIT Parameters class does not allow for names to include special characters.
    This function converts '[', ']', and '.' to their HTML numbers to comply with
    LMFIT.

    Args:
        my_string: A string object

    Returns:
        String object.
    """
    for item in ENCODING_PAIRS:
        my_string = my_string.replace(*item)
    return my_string


def _str_decode(my_string):
    """
    Converts the HTML numbers to '[', ']', and '.' to allow for execution of the
    parameter name to update the simulator.

    Args:
        my_string: A string object

    Returns:
        String Object.
    """
    for item in ENCODING_PAIRS:
        my_string = my_string.replace(*item[::-1])
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
            if key not in EXCLUDE and vals is not None:
                if isinstance(vals, (dict, list)):
                    name_list += _traverse_dictionaries(
                        vals, _str_encode(f"{parent}.{key}")
                    )
                else:
                    name_list += [_str_encode(f"{parent}.{key}")]
    elif isinstance(dictionary, list):
        for i, items in enumerate(dictionary):
            name_list += _traverse_dictionaries(items, _str_encode(f"{parent}[{i}]"))

    return name_list


def _post_sim_LMFIT_params(post_sim):
    """
    Creates an LMFIT Parameters object for SignalProcessor operations
    involved in spectrum fitting

    Args:
        post_sim: SignalProcessor object

    Returns:
        Parameters object
    """
    temp_dict = {}
    # for item in post_sim.operations:
    #     prepend = f"DEP_VAR_{item.dependent_variable}_"
    for i, operation in enumerate(post_sim.operations):
        if isinstance(operation, apo.Gaussian):
            identifier = f"operation_{i}_Gaussian_FWHM"
            arg = operation.FWHM
            temp_dict[f"{identifier}"] = arg
        elif isinstance(operation, apo.Exponential):
            identifier = f"operation_{i}_Exponential_FWHM"
            arg = operation.FWHM
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
    """Updates SignalProcessor operation arguments from an LMFIT Parameters object

    Args:
        params: LMFIT Parameters object
        post_sim: SignalProcessor object
    """
    temp_dict = {}
    arg_dict = {"Gaussian": "FWHM", "Exponential": "FWHM", "Scale": "factor"}
    for param in params:
        # iterating through the parameter list looking for only DEP_VAR
        # (ie post_sim params)
        if "operation_" in param:
            # splitting parameter name to obtain
            # Dependent variable index (var)
            # index of operation in the operation list (opIndex)
            # arg value for the operation (val)
            split_name = param.split("_")
            # var = split_name[split_name.index("VAR") + 1]
            opIndex = split_name[split_name.index("operation") + 1]
            val = params[param].value
            # creating a dictionary of operations and arguments for each dependent
            # variable
            # if f"DepVar_{var}" not in temp_dict.keys():
            #     temp_dict[f"DepVar_{var}"] = {}
            temp_dict[f"{opIndex}_{split_name[2]}"] = val

    # iterate through list of operation lists
    # for item in post_sim.operations:
    # iterating through dictionary with corresponding dependent variable index
    for operation, val in temp_dict.items():
        # creating assignment strings to create the correct address for updating each
        # operation
        split = operation.split("_")
        # dep_var_operation_list = f"post_sim.operations[{item.dependent_variable}]"
        operation_val_update = f"post_sim.operations[{split[0]}].{arg_dict[split[-1]]}"
        assignment = f"={val}"
        exec(operation_val_update + assignment)


def make_LMFIT_parameters(sim, post_sim=None, exclude_key=None):
    """
    Parses the Simulator and PostSimulator objects for a list of LMFIT parameters.
    The parameter name is generated using the following syntax:

    ``sys_i_site_j_attribute1_attribute2``

    for spin system attribute with signature sys[i].sites[j].attribute1.attribute2

    Args:
        sim: a Simulator object.
        post_sim: a SignalProcessor object

    Returns:
        LMFIT Parameters object.
    """
    if not FOUND_LMFIT:
        error = (
            f"The helper function {__name__} requires 'lmfit' module to create lmfit "
            r"paramters. Please install the lmfit module using\n'pip install lmfit'.",
        )
        raise ImportError(error)

    if not isinstance(sim, Simulator):
        raise ValueError(f"Expecting a `Simulator` object, found {type(sim).__name__}.")

    params = Parameters()
    temp_list = _traverse_dictionaries(_list_of_dictionaries(sim.spin_systems))

    # get total abundance scaling factor
    length = len(sim.spin_systems)
    abundance_scale = 100 / sum([sim.spin_systems[i].abundance for i in range(length)])

    # expression for the last abundance.
    last_abund = f"{length - 1}_abundance"
    expression = "-".join([f"{START}{i}_abundance" for i in range(length - 1)])
    expression = "100" if expression == "" else f"100-{expression}"
    for items in temp_list:
        if "_eta" in items:
            params.add(
                name=items, value=eval("sim." + _str_decode(items)), min=0, max=1,
            )
        # last_abund should come before abundance
        elif last_abund in items:
            params.add(
                name=items,
                value=eval("sim." + _str_decode(items)),
                min=0,
                max=100,
                expr=expression,
            )
        elif "abundance" in items:
            params.add(
                name=items,
                value=eval("sim." + _str_decode(items)) * abundance_scale,
                min=0,
                max=100,
            )
        else:
            value = eval("sim." + _str_decode(items))
            params.add(name=items, value=value)

    if post_sim is None:
        return params

    if not isinstance(post_sim, sp.SignalProcessor):
        raise ValueError(
            f"Expecting a `SignalProcessor` object, found {type(post_sim).__name__}."
        )
    temp_params = _post_sim_LMFIT_params(post_sim)
    for item in temp_params:
        params.add(name=item, value=temp_params[item].value)
    # params.add_many(temp_params)

    return params


def LMFIT_min_function(params, sim, post_sim=None):
    """
    The simulation routine to calculate the vector difference between simulation and
    experiment based on the parameters update.

    Args:
        params: Parameters object containing parameters to vary during minimization.
        sim: Simulator object used in the simulation. Initialized with guess fitting
            parameters.
        post_sim: PostSimulator object used in the simulation. Initialized with guess
            fitting parameters.

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
            nameString = "sim." + _str_decode(items)
            executable = f"{nameString} = {values[items]}"
            exec(executable)
        elif "operation_" in items and post_sim is not None:
            _update_post_sim_from_LMFIT_params(params, post_sim)

    sim.run()
    processed_data = post_sim.apply_operations(data=sim.methods[0].simulation)
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
