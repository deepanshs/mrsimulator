# -*- coding: utf-8 -*-
from copy import deepcopy
from typing import ClassVar

from mrsimulator.method import Method
from mrsimulator.method import SpectralDimension
from mrsimulator.utils.error import ImmutableEventError
from mrsimulator.utils.error import NamedMethodError

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


class BaseMethod(Method):
    """BaseMethod class."""

    ndim: ClassVar[int] = 1

    def __init__(self, **kwargs):
        BaseMethod.check(kwargs, self.__class__.ndim)
        super().__init__(**kwargs)

    @classmethod
    def check(cls, kwargs, ndim):
        check_for_number_of_spectral_dimensions(kwargs, ndim)
        if isinstance(kwargs["spectral_dimensions"][0], dict):
            parse_spectral_dimensions(kwargs)
            check_for_atleast_one_events(kwargs)


class Method1D(BaseMethod):
    """Generic one-dimensional spectrum simulation method.

    Example
    -------
    >>> from mrsimulator.methods import Method1D
    >>> method1 = Method1D(
    ...     channels=["87Rb"],
    ...     magnetic_flux_density=7,  # in T
    ...     rotor_angle=54.735 * np.pi / 180,
    ...     rotor_frequency=1e9,
    ...     spectral_dimensions=[
    ...         {
    ...             "count": 1024,
    ...             "spectral_width": 1e4,  # in Hz
    ...             "reference_offset": -4e3,  # in Hz
    ...             "label": "quad only",
    ...             "events": [{"transition_query": [{"P": [-3], "D": [0]}]}],
    ...         }
    ...     ],
    ... )
    """

    ndim: ClassVar[int] = 1
    name: str = "Method1D"
    description: str = "A generic one-dimensional spectrum method."

    def __init__(self, **kwargs):
        # kwargs_copy = deepcopy(kwargs)
        super().__init__(**kwargs)


class Method2D(BaseMethod):
    """Generic two-dimensional spectrum simulation method.

    Example
    -------
    >>> from mrsimulator.methods import Method2D
    >>> method = Method2D(
    ...     channels=["87Rb"],
    ...     magnetic_flux_density=7,  # in T. Global value for `magnetic_flux_density`.
    ...     rotor_angle=0.95531,  # in rads. Global value for the `rotor_angle`.
    ...     spectral_dimensions=[
    ...         {
    ...             "count": 256,
    ...             "spectral_width": 4e3,  # in Hz
    ...             "reference_offset": -5e3,  # in Hz
    ...             "event": [
    ...                 {   # Global value for the `magnetic_flux_density` and
    ...                     # `rotor_angle` is used during this event.
    ...                     "transition_query": {"P": [-3], "D": [0]}
    ...                 }
    ...             ],
    ...         },
    ...         {
    ...             "count": 512,
    ...             "spectral_width": 1e4,  # in Hz
    ...             "reference_offset": -4e3,  # in Hz
    ...             "event": [
    ...                 {   # Global value for `magnetic_flux_density` and user defined
    ...                     # value for `rotor_angle` is used during this event.
    ...                     "rotor_angle": 1.2238,  # in rads
    ...                     "transition_query": {"P": [-1], "D": [0]},
    ...                 }
    ...             ],
    ...         },
    ...     ],
    ...     affine_matrix=[[1, -1], [0, 1]],
    ... )
    """

    ndim: ClassVar[int] = 2
    name: str = "Method2D"
    description: str = "A generic two-dimensional correlation spectrum method."
    rotor_frequency: float = 1.0e12

    def __init__(self, **kwargs):
        kwargs_copy = deepcopy(kwargs)
        # super().check(kwargs_copy, self.__class__.ndim)

        # des = "A generic two-dimensional correlation spectrum method."
        # if "name" not in kwargs_copy:
        #     kwargs_copy["name"] = "Method2D"
        # if "description" not in kwargs_copy:
        #     kwargs_copy["description"] = des
        check_rotor_frequency(kwargs_copy, self.__class__.__name__)
        # if "rotor_frequency" not in kwargs_copy:
        #     kwargs_copy["rotor_frequency"] = 1.0e12

        super().__init__(**kwargs_copy)


class BaseNamedMethod(BaseMethod):
    """BaseNameMethod class."""

    def __init__(self, **kwargs):
        kwargs_copy = deepcopy(kwargs)
        super().check(kwargs_copy, self.__class__.ndim)
        cls_name = self.__class__.__name__
        self.__class__.check_method_compatibility(kwargs_copy)
        check_for_method_name(kwargs_copy, cls_name)
        if self.__class__.ndim > 1:
            check_rotor_frequency(kwargs_copy, cls_name)
        super().__init__(**kwargs_copy)

    @classmethod
    def update(cls, **kwargs):
        return {"spectral_dimensions": [{"events": [{}]} for _ in range(cls.ndim)]}

    @classmethod
    def check_method_compatibility(cls, py_dict):
        """Check for events attribute inside the spectral_dimensions. Events are not
        allowed for NamedMethods."""
        if not isinstance(py_dict["spectral_dimensions"][0], dict):
            return cls.check_when_arg_is_object(py_dict)

        default_method = cls.update(**py_dict)
        default_spectral_dimensions = default_method["spectral_dimensions"]
        for i, item in enumerate(py_dict["spectral_dimensions"]):
            if item["events"] == [{}]:
                item["events"] = default_spectral_dimensions[i]["events"]

            elif item["events"] != default_spectral_dimensions[i]["events"]:
                raise ImmutableEventError(cls.__name__)

        for k, v in default_method.items():
            if k not in py_dict:
                py_dict[k] = v

    @classmethod
    def check_when_arg_is_object(cls, obj_dict):
        default_method = cls.update(**obj_dict)

        py_sp = default_method["spectral_dimensions"]
        obj_sp = obj_dict["spectral_dimensions"]

        for py, obj in zip(py_sp, obj_sp):

            if len(py["events"]) != len(obj.events):
                raise ImmutableEventError(cls.__name__)

            cls.check_event_objects_for_comptibility(py, obj, obj_dict)

    @classmethod
    def check_event_objects_for_comptibility(cls, py, obj, obj_dict):
        required = ["magnetic_flux_density", "rotor_frequency", "rotor_angle"]
        py_obj = SpectralDimension(**py)
        for i, (ev_py, ev_obj) in enumerate(zip(py_obj.events, obj.events)):

            default_obj = SpectralDimension(events=[{}]).events[0]
            obj_keys = ev_obj.dict(exclude={"property_units"}).keys()
            py_keys = py["events"][i].keys()
            for k in obj_keys:
                a = False
                if k in py_keys:
                    a1, a2, a3 = [getattr(_, k) for _ in [ev_obj, default_obj, ev_py]]
                    a = a1 != a2 and a1 != a3 and a2 is not None
                    setattr(ev_obj, k, a3)
                elif k in required and k in obj_dict:
                    a = getattr(ev_obj, k) != obj_dict[k]
                    setattr(ev_obj, k, obj_dict[k])
                if a:
                    raise ImmutableEventError(cls.__name__)


class BaseNamedMethod1D(BaseNamedMethod):
    """Base class for named one-dimensional simulation simulation method."""

    ndim: ClassVar[int] = 1


class BaseNamedMethod2D(BaseNamedMethod):
    """Base class for named one-dimensional simulation simulation method."""

    ndim: ClassVar[int] = 2


class BlochDecaySpectrum(BaseNamedMethod1D):
    """Simulate a Bloch decay spectrum."""

    name: str = "BlochDecaySpectrum"
    description: str = "A one-dimensional Bloch decay spectrum method."

    @classmethod
    def update(cls, **kwargs):
        events = [{"transition_query": [{"ch1": {"P": [-1]}}]}]
        return {
            "spectral_dimensions": [{"events": events}],
        }


class BlochDecayCTSpectrum(BaseNamedMethod1D):
    """Simulate a Bloch decay central transition selective spectrum."""

    name: str = "BlochDecayCTSpectrum"
    description: str = (
        "A one-dimensional central transition selective Bloch decay spectrum method."
    )

    @classmethod
    def update(cls, **kwargs):
        events = [{"transition_query": [{"ch1": {"P": [-1], "D": [0]}}]}]
        return {
            "spectral_dimensions": [{"events": events}],
        }


class BlochDecayCentralTransitionSpectrum(BlochDecayCTSpectrum):
    name: str = "BlochDecayCentralTransitionSpectrum"
    pass


def check_for_number_of_spectral_dimensions(py_dict, n=1):
    """If spectral_dimensions in py_dict, extract and then remove from py_dict."""

    if "spectral_dimensions" not in py_dict:
        py_dict["spectral_dimensions"] = [{} for _ in range(n)]
        return

    m = len(py_dict["spectral_dimensions"])
    if m == n:
        return
    raise ValueError(f"Method requires exactly {n} spectral dimensions, given {m}.")


def parse_spectral_dimensions(py_dict):
    for dim in py_dict["spectral_dimensions"]:
        if "events" in dim.keys():
            for evt in dim["events"]:
                parse_events(evt)
    # return spectral_dimensions


def parse_events(evt):
    if "transition_query" not in evt.keys():
        return

    t_query = evt["transition_query"]
    for item in t_query:
        keys = item.keys()
        if "ch1" not in keys and "ch2" not in keys and "ch3" not in keys:
            item["ch1"] = item


def check_for_atleast_one_events(py_dict):
    for item in py_dict["spectral_dimensions"]:
        if "events" not in item:
            item["events"] = [{}]


def check_for_method_name(py_dict, cls_name):
    if "name" not in py_dict:
        py_dict["name"] = cls_name
        return
    if py_dict["name"] == cls_name:
        return
    name = py_dict["name"]
    raise NamedMethodError(
        f"`name={name} != classname={cls_name}`. Use the class with the same name as "
        "the attribute name or use `Method1D` or `Method2D` class."
    )


def check_rotor_frequency(kwargs, cls_name):
    if "name" in kwargs:
        if kwargs["name"] == "SSB2D":
            return
    if "rotor_frequency" not in kwargs:
        return
    if kwargs["rotor_frequency"] in ["1000000000000.0 Hz", 1.0e12]:
        return
    e = "`rotor_frequency=1e12 Hz` is fixed for 2D Methods and cannot be modified."
    raise ValueError(e)
