# -*- coding: utf-8 -*-
from copy import deepcopy
from os import path

from monty.serialization import loadfn
from mrsimulator.method import Method
from mrsimulator.method.event import Event
from mrsimulator.method.spectral_dimension import SpectralDimension


MODULE_DIR = path.dirname(path.abspath(__file__))
METHODS_DATA = loadfn(path.join(MODULE_DIR, "methods_data.json"))

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


default_units = {
    "magnetic_flux_density": "T",
    "rotor_frequency": "Hz",
    "rotor_angle": "rad",
}


# class namedMethod:
#     def __new__(cls, py_dict):
#         method_name = py_dict["name"]


def prepare_method_structure(template, **kwargs):
    keys = kwargs.keys()
    n_channels = template["number_of_channels"]
    name = template["name"] if "name" not in keys else kwargs["name"]
    desc = (
        template["description"] if "description" not in keys else kwargs["description"]
    )
    label = None if "label" not in keys else kwargs["label"]
    affine_matrix = None if "affine_matrix" not in keys else kwargs["affine_matrix"]

    simulation = None if "simulation" not in keys else kwargs["simulation"]
    experiment = None if "experiment" not in keys else kwargs["experiment"]

    prep = {
        "name": name,
        "description": desc,
        "label": label,
        "simulation": simulation,
        "experiment": experiment,
        "affine_matrix": affine_matrix,
    }

    if "channels" in kwargs:
        prep["channels"] = kwargs["channels"]
        given_n_channels = len(prep["channels"])
        if given_n_channels != n_channels:
            raise ValueError(
                f"The method requires exactly {n_channels} channel(s), "
                f"{given_n_channels} provided."
            )
        kwargs.pop("channels")
    else:
        prep["channels"] = ["1H" for _ in range(n_channels)]
    return prep


def parse_spectral_dimensions(spectral_dimensions, n):
    if spectral_dimensions == [{}]:
        return [{}] * n
    for dim in spectral_dimensions:
        if "events" in dim.keys():
            for evt in dim["events"]:
                parse_events(evt)
    return spectral_dimensions


def parse_events(evt):
    if "transition_query" not in evt.keys():
        return

    t_query = evt["transition_query"]
    if "P" in t_query.keys():
        t_query["P"] = parse_transition_query(t_query["P"])

    if "D" in t_query.keys():
        t_query["D"] = parse_transition_query(t_query["D"])


def parse_transition_query(query):
    if not isinstance(query, list):
        return query
    return {
        "channel-1": [[i] if not isinstance(i, (list, tuple)) else i for i in query]
    }


def generate_method_from_template(template, docstring=""):
    """Generate method object from json template."""
    # constructor
    def constructor(self, spectral_dimensions=[{}], **kwargs):
        # spectral_dimensions_root = deepcopy(spectral_dimensions)
        # kwargs_root = deepcopy(kwargs)
        if isinstance(spectral_dimensions[0], SpectralDimension):
            return Method(spectral_dimensions=spectral_dimensions, **kwargs)

        parse = False
        if "parse" in kwargs:
            parse = kwargs["parse"]
            kwargs.pop("parse")

        local_template = deepcopy(template)
        n_tem = len(local_template["spectral_dimensions"])

        spectral_dimensions = parse_spectral_dimensions(spectral_dimensions, n_tem)
        prep = prepare_method_structure(local_template, **kwargs)
        global_events = local_template["global_event_attributes"]
        ge = set(global_events)
        kw = set(kwargs)
        common = kw.intersection(ge)

        _check_rotor_frequency(common, prep["name"], kwargs)

        n_sp = len(spectral_dimensions)
        if n_tem < n_sp:
            raise ValueError(
                f"The method allows {n_tem} spectral dimension(s), {n_sp} given."
            )
        dim = []
        for i, s in enumerate(local_template["spectral_dimensions"]):
            events = _gen_event(s, spectral_dimensions[i], global_events, parse, kwargs)

            if "events" in spectral_dimensions[i]:
                spectral_dimensions[i].pop("events")

            params = {**spectral_dimensions[i], "events": events}
            params = params if parse else SpectralDimension(**params)
            dim.append(params)

        method = {**prep, "spectral_dimensions": dim}
        method = Method.parse_dict_with_units(method) if parse else Method(**method)
        return method

    @classmethod
    def parse_dict_with_units(cls, py_dict: dict):
        """
        Parse the physical quantities of the method object from as a python dictionary.

        Args:
            py_dict: Dict object.
        """
        return cls.__new__(0, parse=True, **deepcopy(py_dict))

    description = template["description"]
    method = type(
        template["name"],
        (object,),
        {
            "__new__": constructor,
            "__str__": description,
            "__doc__": description + docstring,
            "parse_dict_with_units": parse_dict_with_units,
            "ndim": len(template["spectral_dimensions"]),
        },
    )
    return method


def _fix_strings_in_events(event):
    """Fix the default float values to string with units when parse dict with units
    function is called."""
    for key, val in event.items():
        if key in default_units and isinstance(val, float):
            unit = default_units[key]
            event[key] = f"{val} {unit}"
    return event


def _fill_missing_events_in_template(spectral_dimensions, s_template):
    """Fill the missing events in the template relative to the spectral dimensions."""
    if "events" not in spectral_dimensions:
        return

    s_tem_len = len(s_template["events"])
    sp_evt_len = len(spectral_dimensions["events"])
    if s_tem_len < sp_evt_len:
        diff = sp_evt_len - s_tem_len
        _ = [s_template["events"].append(Event().dict()) for _ in range(diff)]


def _check_rotor_frequency(common, name, kwargs):
    if common == set():
        return

    r_freq = kwargs["rotor_frequency"]
    if "rotor_frequency" in common and r_freq == "1000000000000.0 Hz":
        return

    info = "`, `".join(list(common))
    raise ValueError(f"`{info}` value cannot be modified for {name} method.")


def _gen_event(dim_template, dim_i, global_events, parse, kwargs):

    events = []
    _fill_missing_events_in_template(dim_i, dim_template)

    for j, e in enumerate(dim_template["events"]):
        kw = deepcopy(kwargs)
        ge = deepcopy(global_events)

        # Remove `property_units` from default events instance.
        if "property_units" in e:
            e.pop("property_units")

        if "events" in dim_i:
            kw.update(dim_i["events"][j])

        e.update(kw)
        # prioritize the keyword arguments over the global arguments.
        ge.update(kw)
        e.update(ge)

        params = _fix_strings_in_events(e) if parse else Event(**e)
        events.append(params)

    return events
