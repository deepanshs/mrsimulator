# -*- coding: utf-8 -*-
import mrsimulator.signal_processing as sp
import numpy as np
from lmfit import Parameters
from mrsimulator import Simulator

__author__ = ["Maxwell C Venetos", "Deepansh Srivastava"]
__email__ = ["maxvenetos@gmail.com", "srivastava.89@osu.edu"]

START = "sys_"
ENCODING_PAIRS = [
    ["spin_systems[", START],
    ["].abundance", "_abundance"],
    # site
    ["].sites[", "_site_"],
    ["].isotropic_chemical_shift", "_isotropic_chemical_shift"],
    ["].shielding_symmetric.", "_shielding_symmetric_"],
    ["].quadrupolar.", "_quadrupolar_"],
    # coupling
    ["].couplings[", "_coupling_"],
    ["].isotropic_j", "_isotropic_j"],
    ["].j_symmetric.", "_j_symmetric_"],
    ["].dipolar.", "_dipolar_"],
    # post simulation
    ["].post_simulation", "_POST_SIM_"],
    [".scale", "scale"],
    [".apodization[", "APODIZATION_"],
    ["].args", "_args"],
]

DECODING_PAIRS = [
    ["spin_systems.", "sys_"],
    [".abundance", "_abundance"],
    # site
    [".sites.", "_site_"],
    [".isotropic_chemical_shift", "_isotropic_chemical_shift"],
    [".shielding_symmetric.", "_shielding_symmetric_"],
    [".quadrupolar.", "_quadrupolar_"],
    # coupling
    [".couplings.", "_coupling_"],
    [".isotropic_j", "_isotropic_j"],
    [".j_symmetric.", "_j_symmetric_"],
    [".dipolar.", "_dipolar_"],
]

EXCLUDE = [
    "property_units",
    "isotope",
    "name",
    "label",
    "description",
    "transition_pathways",
]

POST_SIM_DICT = {"Gaussian": "FWHM", "Exponential": "FWHM", "Scale": "factor"}


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
    Parse the string for objects and indexes and return a list.

    Args:
        my_string: A string object

    Returns:
        List of strings with strings resresenting mrsimulator objects and indexes.

    Example:
        >>> string = 'sys_0_site_0_isotropic_chemical_shift'
        >>> _str_decode(string)
        ['spin_systems', '0', 'sites', '0', 'isotropic_chemical_shift']

        >>> string = 'sys_10_abundance'
        >>> _str_decode(string)
        ['spin_systems', '10', 'abundance']
    """
    for item in DECODING_PAIRS:
        my_string = my_string.replace(*item[::-1])
    my_string = my_string.split(".")
    return my_string


def _get_simulator_object_value(sim, string):
    """Parse the string representing the Simulator object dictionary tree format, and
    return its value.

    Args:
        sim: The simulator object.
        string: A string representing the Simulator object dictionary tree format.

    Returns:
        Float. The value of the object.

    Example:
        >>> site = Site(isotropic_chemical_shift=-431)
        >>> sys = SpinSystem(sites=[site], abundance=23)
        >>> sim = Simulator()
        >>> sim.spin_systems.append(sys)

        >>> string = 'sys_0_site_0_isotropic_chemical_shift'
        >>> _get_simulator_object_value(sim, string)
        -431.0

        >>> string = 'sys_0_abundance'
        >>> _get_simulator_object_value(sim, string)
        23.0
    """
    string = _str_decode(string)
    obj = sim
    for attr in string:
        obj = obj[int(attr)] if attr.isnumeric() else obj.__getattribute__(attr)
    return obj


def _update_sim_from_LMFIT_params(sim, string, value):
    """Parse the string representing the Simulator object dictionary tree format, and
    set its value to the input.

    Args:
        sim: The simulator object.
        string: A string representing the Simulator object dictionary tree format.
        value: The value to assign.

    Example:
        >>> site = Site(isotropic_chemical_shift=-431)
        >>> sys = SpinSystem(sites=[site], abundance=23)
        >>> sim = Simulator()
        >>> sim.spin_systems.append(sys)

        >>> string = 'sys_0_site_0_isotropic_chemical_shift'
        >>> _update_sim_from_LMFIT_params(sim, string, 120)
        >>> sim.spin_systems[0].sites[0].isotropic_chemical_shift
        120.0
    """
    string = _str_decode(string)
    obj = sim
    for attr in string[:-1]:
        obj = obj[int(attr)] if attr.isnumeric() else obj.__getattribute__(attr)
    obj.__setattr__(string[-1], value)


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
    if isinstance(dictionary, list):
        return [
            value
            for i, obj in enumerate(dictionary)
            for value in _traverse_dictionaries(obj, _str_encode(f"{parent}[{i}]"))
        ]

    if isinstance(dictionary, dict):
        return [
            item
            for key, value in dictionary.items()
            if key not in EXCLUDE and value is not None
            for item in (
                _traverse_dictionaries(value, _str_encode(f"{parent}.{key}"))
                if isinstance(value, (dict, list))
                else [_str_encode(f"{parent}.{key}")]
            )
        ]

    return []


def _post_sim_LMFIT_params(params, post_sim, index):
    """
    Creates an LMFIT Parameters object for SignalProcessor operations
    involved in spectrum fitting.

    Args:
        post_sim: SignalProcessor object

    Returns:
        Parameters object
    """
    for i, operation in enumerate(post_sim.operations):
        name = operation.__class__.__name__
        if name in POST_SIM_DICT:
            attr = POST_SIM_DICT[name]
            key = f"SP_{index}_operation_{i}_{name}_{attr}"
            val = operation.__getattribute__(attr)
            params.add(name=key, value=val)


def _update_post_sim_from_LMFIT_params(params, post_sim):
    """Updates SignalProcessor operation arguments from an LMFIT Parameters object

    Args:
        params: LMFIT Parameters object
        post_sim: SignalProcessor object
    """
    for param in params:
        # iterating through the parameter list looking for only post_sim params
        if "operation_" in param:
            # splitting parameter name to obtain operations index, operation argument,
            # and its value
            # SP_j_operation_i_function_arg
            split_name = param.split("_")
            sp_idx = int(split_name[1])
            op_idx = int(split_name[3])  # The operation index
            function = split_name[4]
            val = params[param].value  # The value of operation argument parameter

            # update the post_sim object with the parameter updated value.
            post_sim[sp_idx].operations[op_idx].__setattr__(
                POST_SIM_DICT[function], val
            )


def make_LMFIT_parameters(sim: Simulator, post_sim: list = None):
    """An alias of `make_LMFIT_params` function."""
    return make_LMFIT_params(sim, post_sim)


def make_LMFIT_params(sim: Simulator, post_sim: list = None):
    """
    Parses the Simulator and PostSimulator objects for a list of LMFIT parameters.
    The parameter name is generated using the following syntax:

    ``sys_i_site_j_attribute1_attribute2``

    for spin system attribute with signature sys[i].sites[j].attribute1.attribute2

    Args:
        sim: A Simulator object.
        post_sim: A list of SignalProcessor object of length equal to the length of
            methods in the Simulator object.

    Returns:
        LMFIT Parameters object.
    """
    params = Parameters()
    make_simulator_params(params, sim)

    if post_sim is not None:
        params.update(make_signal_processor_params(post_sim))

    return params


def make_signal_processor_params(post_sim):
    post_sim = post_sim if isinstance(post_sim, list) else [post_sim]

    params = Parameters()
    for i, processor in enumerate(post_sim):
        if not isinstance(processor, sp.SignalProcessor):
            name = type(processor).__name__
            raise ValueError(f"Expecting a `SignalProcessor` object, found {name}.")
        _post_sim_LMFIT_params(params, processor, i)
    return params


def make_simulator_params(params, sim):
    if not isinstance(sim, Simulator):
        raise ValueError(f"Expecting a `Simulator` object, found {type(sim).__name__}.")

    temp_list = _traverse_dictionaries(_list_of_dictionaries(sim.spin_systems))

    # get total abundance scaling factor
    length = len(sim.spin_systems)
    abundance_scale = 100 / sum([sim.spin_systems[i].abundance for i in range(length)])

    # expression for the last abundance.
    last_abundance = f"{length - 1}_abundance"
    expression = "-".join([f"{START}{i}_abundance" for i in range(length - 1)])
    expression = "100" if expression == "" else f"100-{expression}"
    for items in temp_list:
        value = _get_simulator_object_value(sim, items)
        if "_eta" in items:
            params.add(name=items, value=value, min=0, max=1)

        # last_abundance should come before abundance
        elif last_abundance in items:
            params.add(name=items, value=value, min=0, max=100, expr=expression)

        elif "abundance" in items:
            params.add(name=items, value=value * abundance_scale, min=0, max=100)
        else:
            params.add(name=items, value=value)


def LMFIT_min_function(
    params: Parameters, sim: Simulator, post_sim: list = None, sigma: list = None
):
    """
    The simulation routine to calculate the vector difference between simulation and
    experiment based on the parameters update.

    Args:
        params: Parameters object containing parameters for OLS minimization.
        sim: Simulator object.
        post_sim: A list of PostSimulator objects corresponding to the methods in the
            Simulator object.
        sigma: A list of standard deviations corresponding to the experiments in the
            Simulator.methods attribute
    Returns:
        Array of the differences between the simulation and the experimental data.
    """
    post_sim = post_sim if isinstance(post_sim, list) else [post_sim]
    if sigma is None:
        sigma = [1.0 for _ in sim.methods]
    sigma = sigma if isinstance(sigma, list) else [sigma]

    values = params.valuesdict()
    for items in values:
        if "operation_" not in items:
            _update_sim_from_LMFIT_params(sim, items, values[items])
        elif "operation_" in items and post_sim is not None:
            _update_post_sim_from_LMFIT_params(params, post_sim)

    sim.run()

    processed_data = [
        item.apply_operations(data=data.simulation)
        for item, data in zip(post_sim, sim.methods)
    ]

    diff = np.asarray([])
    for processed_datum, mth, sigma_ in zip(processed_data, sim.methods, sigma):
        datum = 0
        for decomposed_datum in processed_datum.y:
            datum += decomposed_datum.components[0].real

        # If data has negative increment, reverse the data before taking the difference.
        exp_data = mth.experiment.y[0].components[0]
        index = [
            -i - 1 for i, x in enumerate(mth.experiment.x) if x.increment.value < 0
        ]
        exp_data = exp_data if index == [] else np.flip(exp_data, axis=tuple(index))
        diff = np.append(diff, (exp_data - datum) / sigma_)

    return diff
