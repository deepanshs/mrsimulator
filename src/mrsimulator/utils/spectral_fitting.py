# -*- coding: utf-8 -*-
import mrsimulator.signal_processing as sp
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
    [".shielding_j.", "_shielding_j_"],
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


def _set_simulator_object_value(sim, string, value):
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
        >>> _set_simulator_object_value(sim, string, 120)
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
    involved in spectrum fitting.

    Args:
        post_sim: SignalProcessor object

    Returns:
        Parameters object
    """
    params = Parameters()

    for i, operation in enumerate(post_sim.operations):
        name = operation.__class__.__name__
        if name in POST_SIM_DICT:
            attr = POST_SIM_DICT[name]
            key = f"operation_{i}_{name}_{attr}"
            val = operation.__getattribute__(attr)
            params.add(name=key, value=val)

    return params


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
            split_name = param.split("_")
            opIndex = int(split_name[1])  # The operation index
            val = params[param].value  # The value of operation argument parameter

            # update the post_sim object with the parameter updated value.
            post_sim.operations[opIndex].__setattr__(POST_SIM_DICT[split_name[2]], val)


def make_LMFIT_parameters(sim, post_sim=None, exclude_key=None):
    """An alias of `make_LMFIT_params` function."""
    return make_LMFIT_params(sim, post_sim, exclude_key)


def make_LMFIT_params(sim, post_sim=None, exclude_key=None):
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
    if not isinstance(sim, Simulator):
        raise ValueError(f"Expecting a `Simulator` object, found {type(sim).__name__}.")

    params = Parameters()
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


def LMFIT_min_function(params, sim, post_sim=None, sigma=1):
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
    # if not isinstance(params, Parameters):
    #     raise ValueError(
    #         f"Expecting a `Parameters` object, found {type(params).__name__}."
    #     )
    # if not isinstance(post_sim, sp.SignalProcessor) or post_sim is None:
    #     raise ValueError(
    #         f"Expecting a `SignalProcessor` object, found {type(post_sim).__name__}."
    #     )

    # if not isinstance(sim, Simulator):
    #  raise ValueError(f"Expecting a `Simulator` object, found {type(sim).__name__}.")

    values = params.valuesdict()
    for items in values:
        if "operation_" not in items:
            _set_simulator_object_value(sim, items, values[items])
        elif "operation_" in items and post_sim is not None:
            _update_post_sim_from_LMFIT_params(params, post_sim)

    sim.run()
    processed_data = post_sim.apply_operations(data=sim.methods[0].simulation)

    datum = 0
    if sim.config.decompose_spectrum == "spin_system":
        for decomposed_datum in processed_data.y:
            datum += decomposed_datum.components[0].real
    else:
        datum = processed_data.y[0].components[0].real

    diff = sim.methods[0].experiment.y[0].components[0] - datum
    return diff / sigma

    # MULTIPLE EXPERIMENTS
    # for i, method in enumerate(sim.methods):
    #     y_factored = method.apodize().real
    #     residual = np.append(
    #         residual,
    #         method.experiment.dependent_variables[0].components[0].real - y_factored,
    #     )

    # return residual
