from copy import deepcopy
from typing import ClassVar

from mrsimulator.method import Method
from mrsimulator.method import SpectralDimension
from mrsimulator.utils.error import ImmutableEventError
from mrsimulator.utils.error import NamedMethodError
from pydantic import Field
from pydantic import PrivateAttr
from pydantic import validator

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


class BaseNamedMethod(Method):
    """BaseNamedMethod class."""

    _named_method: bool = PrivateAttr(True)
    _num_channels: ClassVar[int] = 1
    ndim: ClassVar[int] = 1

    def __init__(self, **kwargs):
        kwargs_copy = deepcopy(kwargs)
        Method.check(kwargs_copy, is_named_method=True, ndim=self.__class__.ndim)
        self.__class__.check_method_compatibility(kwargs_copy)
        super().__init__(**kwargs_copy)

    @validator("name", pre=True, always=True)
    def check_for_method_name(cls, v, *, values, **kwargs):
        if v == cls.__name__:
            return v
        raise NamedMethodError(v, cls.__name__)

    @validator("rotor_frequency", pre=True, always=True)
    def check_rotor_frequency(cls, v, *, values, **kwargs):
        a = (
            True if "name" not in values else values["name"] == "SSB2D",
            v in ["1000000000000.0 Hz", 1.0e12],
            cls.ndim == 1,
        )
        if any(a):
            return v
        raise ValueError(
            "`rotor_frequency=1e12 Hz` is fixed for all 2D named Methods, except SSB2D,"
            " and cannot be modified."
        )

    @validator("channels", pre=True, always=True)
    def check_channel_count(cls, v, *, values, **kwargs):
        # NOTE: This check will need to change if multi-isotope named methods added
        if len(v) != cls._num_channels:
            raise ValueError(
                f"{cls.__name__} only supports {cls._num_channels} channel(s). "
                f"Got {len(v)} channels."
            )
        return v

    @classmethod
    def update(cls, **kwargs):
        return {"spectral_dimensions": [{"events": [{}]} for _ in range(cls.ndim)]}

    @classmethod
    def check_method_compatibility(cls, py_dict):
        """Check for events attribute inside the spectral_dimensions. Events are not
        allowed for NamedMethods."""
        sp_list = py_dict["spectral_dimensions"]
        check_SD = [isinstance(sp, SpectralDimension) for sp in sp_list]
        if all(check_SD):
            return cls.check_when_arg_is_object(py_dict)

        default_method = cls.update(**py_dict)
        default_sp_list = default_method["spectral_dimensions"]

        for i, item in enumerate(sp_list):
            # If no events in SpectralDimension, set to default events
            if "events" not in item or item["events"] == [{}] or item["events"] == []:
                item["events"] = default_sp_list[i]["events"]

            elif item["events"] != default_sp_list[i]["events"]:
                raise ImmutableEventError(cls.__name__)

        for k, v in default_method.items():
            if k not in py_dict:
                py_dict[k] = v

    @classmethod
    def check_when_arg_is_object(cls, method_dict):
        default_method = cls.update(**method_dict)
        default_sp = default_method["spectral_dimensions"]
        obj_sp = method_dict["spectral_dimensions"]

        for i, (dflt_dim, obj_dim) in enumerate(zip(default_sp, obj_sp)):
            if len(dflt_dim["events"]) != len(obj_dim.events) and obj_dim.events != []:
                raise ImmutableEventError(cls.__name__)

            if obj_dim.events == []:
                obj_sp[i] = obj_dim.json(units=False)
                if "events" not in obj_sp[i]:
                    obj_sp[i]["events"] = dflt_dim["events"]
                obj_sp[i] = SpectralDimension(**obj_sp[i])

            cls.check_event_objects_for_compatibility(dflt_dim, obj_dim, method_dict)

        for k, v in default_method.items():
            if k not in method_dict:
                method_dict[k] = v

    @classmethod
    def check_event_objects_for_compatibility(cls, default_dim, obj_dim, method_dict):
        """Checks Events for compatibility and sets global method attributes

        Args:
            default_dim (dict): Dict representation of SpectralDimension in base method
            obj_dim (SpectralDimension): User-passed SpectralDimension object to check
            method_dict (dict): Dict representation of passed method
        """
        required = ["magnetic_flux_density", "rotor_frequency", "rotor_angle"]
        check_dim = SpectralDimension(**default_dim)
        for i, (ev_check, ev_obj) in enumerate(zip(check_dim.events, obj_dim.events)):

            default_obj = SpectralDimension(events=[{}]).events[0]
            obj_keys = ev_obj.dict(exclude={"property_units"}).keys()
            check_keys = default_dim["events"][i].keys()
            for k in obj_keys:  # iterate over event attributes
                fail = False
                if k in check_keys:
                    obj_attr, default_attr, check_attr = (
                        getattr(_, k) for _ in [ev_obj, default_obj, ev_check]
                    )
                    fail_1 = obj_attr != default_attr  # not default (user passed value)
                    fail_2 = obj_attr != check_attr  # passed attr does not match base
                    fail_3 = default_attr is not None
                    fail = fail_1 and fail_2 and fail_3
                    setattr(ev_obj, k, check_attr)
                elif k in required and k in method_dict:
                    # True if passed attr does not match global attr defined by method
                    fail = getattr(ev_obj, k) != method_dict[k]
                    # Set event attr to global method attr
                    setattr(ev_obj, k, method_dict[k])
                if fail:
                    raise ImmutableEventError(cls.__name__)


class BaseNamedMethod1D(BaseNamedMethod):
    """Base class for named one-dimensional simulation simulation method."""

    ndim: ClassVar[int] = 1


class BaseNamedMethod2D(BaseNamedMethod):
    """Base class for named two-dimensional simulation simulation method."""

    ndim: ClassVar[int] = 2
    rotor_frequency: float = Field(default=1.0e12, ge=0.0)


class BlochDecaySpectrum(BaseNamedMethod1D):
    """Simulate a Bloch decay spectrum."""

    name: str = "BlochDecaySpectrum"
    description: str = "A one-dimensional Bloch decay spectrum method."

    class Config:
        extra = "forbid"

    @classmethod
    def update(cls, **kwargs):
        events = [{"transition_queries": [{"ch1": {"P": [-1]}}]}]
        return {
            "spectral_dimensions": [{"events": events}],
        }


class BlochDecayCTSpectrum(BaseNamedMethod1D):
    """Simulate a Bloch decay central transition selective spectrum."""

    name: str = "BlochDecayCTSpectrum"
    description: str = (
        "A one-dimensional central transition selective Bloch decay spectrum method."
    )

    class Config:
        extra = "forbid"

    @classmethod
    def update(cls, **kwargs):
        events = [{"transition_queries": [{"ch1": {"P": [-1], "D": [0]}}]}]
        return {
            "spectral_dimensions": [{"events": events}],
        }


# Class Aliases
class BlochDecayCentralTransitionSpectrum(BlochDecayCTSpectrum):
    name: str = "BlochDecayCentralTransitionSpectrum"

    class Config:
        extra = "forbid"

    def __init__(self, **kwargs):
        DeprecationWarning(
            "BlochDecayCentralTransitionSpectrum is deprecated, use ",
            "BlochDecayCTSpectrum class instead",
        )
        super().__init__(**kwargs)
