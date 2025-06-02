import csdmpy as cp
import mrsimulator.signal_processor as sp
import numpy as np
from lmfit import Parameters
from mrsimulator import Simulator
from mrsimulator.models.utils import LineShapeKernel

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

POST_SIM_DICT = {
    "Gaussian": {"FWHM": "FWHM"},
    "Exponential": {"FWHM": "FWHM"},
    "Scale": {"factor": "factor"},
    "ConstantOffset": {"offset": "offset"},
    "Linear": {"amplitude": "amplitude", "offset": "offset"},
}


def _str_encode(my_string):
    """LMFIT Parameters class does not allow for names to include special characters.
    This function replaces '[', ']', and '.' to '_' to comply with LMFIT rules.

    Args:
        my_string: A string object

    Returns:
        String object.
    """
    for item in ENCODING_PAIRS:
        my_string = my_string.replace(*item)
    return my_string


def _str_decode(my_string):
    """Parse the string for objects and indexes and return a list.

    Args:
        my_string: A string object

    Returns:
        List of strings with strings representing mrsimulator objects and indexes.

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


def _list_of_dictionaries(my_list):
    """Helper function for traverse_dictionaries function which will return a list of
    dictionaries.

    Args:
        my_list: A list object

    Returns:
        List Object.
    """
    return [item.dict() for item in my_list]


def _traverse_dictionaries(instance, parent="spin_systems"):
    """Parses through the instance object contained within the parent object and return
    a list of attributes that are populated.

    Args:
        instance: An instance object from the parent object.
        parent: a string object used to create the addresses of the SpinSystem
            attributes.

    Returns:
        List Object.
    """
    if isinstance(instance, list):
        return [
            value
            for i, obj in enumerate(instance)
            for value in _traverse_dictionaries(obj, _str_encode(f"{parent}[{i}]"))
        ]

    if isinstance(instance, dict):
        return [
            item
            for key, value in instance.items()
            if key not in EXCLUDE and value is not None
            for item in (
                _traverse_dictionaries(value, _str_encode(f"{parent}.{key}"))
                if isinstance(value, (dict, list))
                else [_str_encode(f"{parent}.{key}")]
            )
        ]

    return []


def _make_params_single_processor(params, process, index):
    """Creates an LMFIT Parameters object for SignalProcessor operations involved in
    spectrum fitting.

    Args:
        params: LMFIT parameters object.
        process: SignalProcessor object at index *index*.
        int index: List index of the SignalProcessor object.

    Returns:
        Parameters object.
    """
    _ = [
        params.add(
            name=f"SP_{index}_operation_{i}_{operation.__class__.__name__}_{attr}",
            value=getattr(operation, attr),
        )
        for i, operation in enumerate(process.operations)
        if operation.__class__.__name__ in POST_SIM_DICT
        for attr in POST_SIM_DICT[operation.__class__.__name__]
    ]


def make_signal_processor_params(processors: list):
    """Parse the list of SignalProcessor objects for a list of LMFIT parameters.

    Args:
        processors: List of SignalProcessor objects. The order of the list must match
            the order of the methods in the Simulator object.
    """
    processors = processors if isinstance(processors, list) else [processors]

    correct_type = all([isinstance(obj, sp.SignalProcessor) for obj in processors])
    if not correct_type:
        raise ValueError("Expecting a list of `SignalProcessor` objects.")

    params = Parameters()
    _ = [
        _make_params_single_processor(params, proc, i)
        for i, proc in enumerate(processors)
    ]
    return params


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
        obj = obj[int(attr)] if attr.isnumeric() else getattr(obj, attr)
    return obj


def make_simulator_params(
    sim: Simulator = Simulator(), spin_system_models: list = [], include={}
):
    """Parse the Simulator object for a list of LMFIT parameters.

    Args:
        Simulator sim: Simulator object.
    """
    params = Parameters()
    if not isinstance(sim, Simulator):
        raise ValueError(f"Expecting a `Simulator` object, found {type(sim).__name__}.")

    temp_list = _traverse_dictionaries(_list_of_dictionaries(sim.spin_systems))

    # get total abundance scaling factor
    sys_length = len(sim.spin_systems)
    sys_model_length = len(spin_system_models)
    total_abundance = 0.0
    total_abundance += sum(sys.abundance for sys in sim.spin_systems)
    total_abundance += sum(sys.abundance for sys in spin_system_models)
    abundance_scale = 100.0 / total_abundance

    # expression for the last abundance.
    if sys_length > 0:
        last_abundance = f"{sys_length - 1}_abundance"
        expression = "-".join([f"{START}{i}_abundance" for i in range(sys_length - 1)])
        expression = "100" if expression == "" else f"100-{expression}"

    skip_last = sys_length > 0
    if sys_model_length > 0:
        param_dist = make_distribution_params(
            spin_system_models, norm=abundance_scale, skip_last=skip_last
        )
        _ = params.update(param_dist)

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

    get_simulator_method_parameters(params, sim, include)
    return params


def get_simulator_method_parameters(params: Parameters, sim: Simulator, include={}):
    """LMFIT parameters for the method attribute in `include`

    Args:
        Parameters params: LMFIT parameters object.
        Simulator sim: Simulator object.
        set include: Set of method attributes to include.
    """
    if "rotor_frequency" in include:
        values = [
            [
                ev.rotor_frequency
                for sp in mth.spectral_dimensions
                for ev in sp.events
                if ev.rotor_frequency != 1e12
            ]
            for i, mth in enumerate(sim.methods)
            if mth._named_method
        ]

        _ = [
            params.add(
                name=f"mth_{i}_rotor_frequency",
                value=val[0],
                min=val[0] - 100,
                max=val[0] + 100,
            )
            for i, val in enumerate(values)
            if val != []
        ]

    return params


def make_LMFIT_parameters(sim: Simulator, processors: list = None, include={}):
    """An alias of `make_LMFIT_params` function."""
    return make_LMFIT_params(sim, processors)


def make_LMFIT_params(
    sim: Simulator = Simulator(),
    processors: list = None,
    spin_system_models: list = [],
    include={},
):
    r"""Parse the Simulator and PostSimulator objects for a list of LMFIT parameters.

    Args:
        Simulator sim: Simulator object.
        list processors: List of SignalProcessor objects. The order must match the order
            of methods within the simulator object.
        set include: Set of keywords from the method object to include as a fitting
            parameter. Default is {}.

    The parameter name associated with the spin system within the Simulator object is
    generated using the following nomenclature- *sys_i_site_j_attribute1_attribute2*
    for attribute with signature ``sim.spin_systems[i].sites[j].attribute1.attribute2``

    Here, *sys_i* refers to the spin system at index *i*, *site_j* refers to the site at
    index *j* within the :math:i^\text{th} spin system, and *attribute1* and
    *attribute2* are the site attributes.

    **For examples:**

    ``sim.spin_systems[1].sites[0].isotropic_chemical_shift`` parametrizes to
    *sys_1_site_0_isotropic_chemical_shift* while
    ``sim.spin_systems[0].sites[1].quadrupolar.Cq`` to *sys_0_site_1_quadrupolar_Cq*.

    Returns:
        LMFIT Parameters object.
    """
    params = Parameters()
    params.update(
        make_simulator_params(
            sim=sim, spin_system_models=spin_system_models, include=include
        )
    )

    proc = make_signal_processor_params(processors) if processors is not None else None
    _ = params.update(proc) if proc is not None else None

    return params


def _update_simulator_from_LMFIT_params(params, sim: Simulator):
    """Parse the string representing the Simulator object dictionary tree format, and
    set its value to the input.

    Args:
        sim: The simulator object.
        string: A string representing the Simulator object dictionary tree format.
        value: The value to assign.

    Example:
        if params['sys_i_site_j_isotropic_chemical_shift'].value = 120
        the function sets, sim.spin_systems[i].sites[j].isotropic_chemical_shift = 120.0
    """

    def set_sys_value(obj, key, value):
        ids = _str_decode(key)
        for attr in ids[:-1]:
            obj = obj[int(attr)] if attr.isnumeric() else obj.__dict__[attr]
        obj.__dict__[ids[-1]] = value

    def set_mth_value(obj, key, value):
        index = int(key.split("_")[1])
        _ = [
            setattr(sp.events[0], "rotor_frequency", value)
            for sp in obj.__dict__["methods"][index].spectral_dimensions
            if sp.events[0].rotor_frequency != 1e12
        ]

    values = params.valuesdict()
    _ = [set_sys_value(sim, k, v) for k, v in values.items() if "sys_" in k]
    _ = [set_mth_value(sim, k, v) for k, v in values.items() if "mth_" in k]


def _update_processors_from_LMFIT_params(params, processors: list):
    """Updates SignalProcessor operation arguments from an LMFIT Parameters object.

    Args:
        params: LMFIT Parameters object.
        processors: list of SignalProcessor objects.

    Example:
        if params['SP_i_operations_j_Scale_factor'].value = 10
        the function sets, processors[i].operations[j].Scale.factor = 10.0
    """

    def set_value(obj, key, value):
        ids = key.split("_")
        # The signal processor index, operation indexes, and function argument
        sp, op, arg = int(ids[1]), int(ids[3]), POST_SIM_DICT[ids[4]][ids[5]]
        obj[sp].__dict__["operations"][op].__dict__[arg] = value

    values = params.valuesdict()
    _ = [set_value(processors, k, v) for k, v in values.items() if "operation_" in k]


def update_mrsim_obj_from_params(params, sim: Simulator, processors: list = None):
    """Update the mrsimulator Simulator and SignalProcessor objects from the LMFIT
    Parameters obj

    Args:
        params: Parameters object containing parameters for OLS minimization.
        sim: Simulator object.
        processors: A list of SignalProcessor objects corresponding to the methods in
            the Simulator object.
    """
    _update_simulator_from_LMFIT_params(params, sim)
    _update_processors_from_LMFIT_params(params, processors)


# def get_correct_data_order(method):
#     """If data has negative increment, reverse the data."""
#     if "experiment" in method._metadata:
#         return method._metadata["experiment"]
#     data = method.experiment
#     y = data.y[0].components[0]
#     idx = [-i - 1 for i, x in enumerate(data.x) if x.increment.value < 0]
#     method._metadata["experiment"] = y if idx == [] else np.flip(y, axis=tuple(idx))
#     return method._metadata["experiment"]


def _check_for_experiment_data(methods_list: list):
    """Ensures all Method objects in the passed list contain experimental data,
    otherwise, a ValueError is raised.

    Args:
        (list) methods_list: A list of Method objects

    Raises:
        ValueError if a method does not contain experimental data
    """
    # Get index of Method only if Method does not hold experimental data
    no_exp = [i for i, mth in enumerate(methods_list) if mth.experiment is None]

    if no_exp:  # At least one item in the list
        raise ValueError(
            f"No experimental data found for the method at index {no_exp}. "
            "Attach an experimental dataset to each method before fitting."
        )


def LMFIT_min_function(
    params: Parameters,
    sim: Simulator,
    processors: list = None,
    sigma: list = None,
    **sim_kwargs,
):
    """The simulation routine to calculate the vector difference between simulation and
    experiment based on the parameters update.

    Args:
        params: Parameters object containing parameters for OLS minimization.
        sim: Simulator object.
        processors: A list of :py:class:~`mrsimulator.signal_processor.Processor`
            objects to apply post-simulation processing to the simulated spectra.
        sigma: A list of standard deviations corresponding to the experiments in the
            :py:attr:~`mrsimulator.Simulator.methods` attribute.
        sim_kwargs: Keyword arguments to pass to the
            :py:mth:~`mrsimulator.Simulator.run()` method.
    Returns:
        Array of the differences between the simulation and the experimental datasets.
    """
    processors = processors if isinstance(processors, list) else [processors]
    sigma = [1.0 for _ in sim.methods] if sigma is None else sigma
    sigma = sigma if isinstance(sigma, list) else [sigma]
    sigma = [item.value if hasattr(item, "value") else item for item in sigma]

    _check_for_experiment_data(sim.methods)
    update_mrsim_obj_from_params(params, sim, processors)

    sim.run(**sim_kwargs)

    processed_dataset = [
        item.apply_operations(dataset=data.simulation)
        for item, data in zip(processors, sim.methods)
    ]

    diff = np.asarray([])
    for processed_datum, mth, sigma_ in zip(processed_dataset, sim.methods, sigma):
        datum = 0
        for decomposed_datum in processed_datum.y:
            datum += decomposed_datum.components[0].real

        exp_data = mth.experiment.y[0].components[0]
        diff = np.append(diff, (exp_data - datum) / sigma_)
    return diff


def bestfit(sim: Simulator, processors: list = None):
    """Return a list of best-fit spectrum ordered relative to the methods in the
    simulator object.

    Args:
        Simulator sim: The simulator object.
        list processors: List of SignalProcessor objects ordered according to the
            methods in the simulator object.
    """
    processors = processors if isinstance(processors, list) else [processors]
    sim.run()

    return [
        proc.apply_operations(dataset=mth.simulation).real
        for mth, proc in zip(sim.methods, processors)
    ]


def add_csdm_dvs(data):
    new_data = data.split()
    new_csdm = 0
    for item in new_data:
        new_csdm += item
    return new_csdm if new_data != [] else None


def residuals(sim: Simulator, processors: list = None):
    """Return a list of residuals corresponding to the best-fit spectrum. The list is
    based on the order of methods in the simulator object.

    Args:
        Simulator sim: The simulator object.
        list processors: List of SignalProcessor objects ordered according to the
            methods in the simulator object."""
    fits = bestfit(sim, processors)
    residual_ = [add_csdm_dvs(item) for item in fits]

    for res, mth in zip(residual_, sim.methods):
        exp_data = mth.experiment.y[0].components[0]
        res.y[0].components[0] -= exp_data
        res.y[0].components[0] *= -1

    return residual_


def _apply_iso_shift(csdm_obj, iso_shift_ppm, larmor_freq):
    """Apply isotropic chemical shift to a CSDM object using the FFT shift theorem."""
    csdm_obj = csdm_obj.fft()
    time_coords = csdm_obj.x[0].coordinates.to("s").value
    iso_shift_Hz = larmor_freq * iso_shift_ppm
    csdm_obj.y[0].components[0] *= np.exp(-np.pi * 2j * time_coords * iso_shift_Hz)
    csdm_obj = csdm_obj.fft()

    return csdm_obj


def make_distribution_params(
    spin_system_models: list, norm: float, skip_last: bool = False
) -> Parameters:
    """Generate LMfit Parameters object for spin system distribution model.
    The distribution has the following parameters:
    """
    n_dist = len(spin_system_models)

    expr_terms = []
    params = Parameters()
    for i, model in enumerate(spin_system_models):
        # normalize the abundance
        model.abundance *= norm
        params = model.add_lmfit_params(params, i)

        # Set last abundance parameter as a non-free variable
        model_prefix = model.model_name
        if i < n_dist - 1 or skip_last:
            expr_terms.append(f"{model_prefix}_{i}_abundance")
        else:
            expr = "-".join(expr_terms)
            expr = "100" if expr == "" else f"100 - {expr}"
            params[f"{model_prefix}_{i}_abundance"].vary = False
            params[f"{model_prefix}_{i}_abundance"].expr = expr
    return params


def _generate_distribution_spectrum(
    params: Parameters,
    kernel: LineShapeKernel,
    spin_system_models: list,
    processor: sp.SignalProcessor = None,
) -> cp.CSDM:
    """Helper function for generating a spectrum from a set of LMfit Parameters and
    and arguments for defining the grid, kernel, etc. The other functions used in the
    least-squares minimization uses this function to reduce code overlap.

    Arguments:
        (Parameters) params: The LMfit parameters object holding parameters used during
            the least-squares minimization.
        (np.ndarray) kernel The pre-computed lineshape kernel.
        (list) spin_system_models:  A list of spin system model.
        (sp.SignalProcessor) processor: A
            :py:class:~`mrsimulator.signal_processor.Processor` object used to apply
            post-simulation signal processing to the resulting spectrum.
            TODO: expand processor to apply to individual distributions (as a list)

    Returns:
        Guess spectrum as a cp.CSDM object

    """
    method = kernel.method
    isotope = method.channels[0]
    larmor_freq = isotope.larmor_freq(B0=method.magnetic_flux_density) / 1.0e6  # MHz
    larmor_freq *= isotope.ref_larmor_ratio
    exp_spectrum = method.experiment

    guess_spectrum = exp_spectrum.copy()
    guess_spectrum.y[0].components[:] = 0  # Initialize guess spectrum with zeros
    dv = cp.as_dependent_variable(np.empty(guess_spectrum.y[0].components.size))

    for i, model in enumerate(spin_system_models):
        model.update_lmfit_params(params, i)
        _, _, amp = model.pdf(kernel.pos)

        # Dot amplitude with kernel, then package as CSDM object
        spec_tmp = cp.CSDM(dimensions=exp_spectrum.x, dependent_variables=[dv])
        spec_tmp.y[0].components[0] = np.dot(amp.ravel(), kernel.kernel)

        # Apply isotropic chemical shift to distribution using FFT shift theorem
        spec_tmp = _apply_iso_shift(
            csdm_obj=spec_tmp,
            iso_shift_ppm=model.mean_isotropic_chemical_shift,
            larmor_freq=larmor_freq,
        ).real
        spec_tmp *= model.abundance
        guess_spectrum.y[0].components[0] += spec_tmp.y[0].components[0]

    if processor is not None:
        _update_processors_from_LMFIT_params(params, [processor])
        guess_spectrum = processor.apply_operations(dataset=guess_spectrum).real

    return guess_spectrum


def LMFIT_min_function_dist(
    params: Parameters,
    kernel: LineShapeKernel,
    spin_system_models: list,
    processor: sp.SignalProcessor = None,
) -> np.ndarray:
    """The minimization routine for fitting a set of Czjzek models to an experimental
    spectrum.

    Arguments:
        (Parameters) params: The LMfit parameters object holding parameters used during
            the least-squares minimization.
        (np.ndarray) kernel: The pre-computed lineshape kernel.
        (list) spin_system_models: A list of spin system distribution models.
        (sp.SignalProcessor) processor: A
            :py:class:~`mrsimulator.signal_processor.Processor` object used to apply
            post-simulation signal processing to the resulting spectrum.
            TODO: expand processor to apply to individual distributions (as a list)

        Returns:
            A residual array as a numpy array.
    """
    guess_spectrum = _generate_distribution_spectrum(
        params, kernel, spin_system_models, processor
    )
    exp_spectrum = kernel.method.experiment
    return (exp_spectrum - guess_spectrum).y[0].components[0]


def bestfit_dist(
    params: Parameters,
    kernel: LineShapeKernel,
    spin_system_models: list,
    processor: sp.SignalProcessor = None,
) -> cp.CSDM:
    """Returns the best-fit spectrum as a CSDM object"""
    return _generate_distribution_spectrum(
        params, kernel, spin_system_models, processor
    )


def residuals_dist(
    params: Parameters,
    kernel: LineShapeKernel,
    spin_system_models: list,
    processor: sp.SignalProcessor = None,
) -> cp.CSDM:
    """Returns the residuals spectrum as a CSDM object"""
    bestfit_ = _generate_distribution_spectrum(
        params, kernel, spin_system_models, processor
    )

    exp_spectrum = kernel.method.experiment
    return exp_spectrum - bestfit_
