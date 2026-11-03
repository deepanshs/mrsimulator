import numpy as np
from astropy.units import Quantity
from csdmpy.units import string_to_quantity


__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"

CONST = string_to_quantity("1")


def _get_broadcast_shape(array, dim, ndim):
    """Return the broadcast shape of a vector `array` at dimension `dim` for `ndim`
    total dimensions."""
    none = [np.newaxis for _ in range(ndim + 1)]
    dim = [dim] if isinstance(dim, int) else dim
    for dim_ in dim:
        none[-dim_ - 1] = slice(None, None, None)
    return array[tuple(none)]


def _str_to_quantity(v, values, prop_name):
    if isinstance(v, str):
        quantity = string_to_quantity(v)
        # keep the units of other properties set from strings, but not the
        # dimensionless defaults (CONST)
        property_units = {
            k: u
            for k, u in (values.get("property_units", None) or {}).items()
            if not isinstance(u, Quantity)
        }
        property_units[prop_name] = quantity.unit
        values["property_units"] = property_units
        return quantity.value
    if isinstance(v, float):
        return v
