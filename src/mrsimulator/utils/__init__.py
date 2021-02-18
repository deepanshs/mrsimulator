# -*- coding: utf-8 -*-

__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"


def get_spectral_dimensions(csdm_object):
    """
    Extract the count, spectral_width, and reference_offset parameters, associated with
    the spectral dimensions of the method, from the CSDM dimension objects.

    Args:
        csdm_object: A CSDM object holding the measurement dataset.

    Returns:
        A list of dict objects, where each dict containts the count, spectral_width, and
        reference_offset.
    """
    result = []
    for dim in csdm_object.dimensions:
        count = dim.count
        increment = dim.increment.to("Hz").value
        ref = dim.coordinates_offset.to("Hz").value
        sw = count * increment
        co = ref if dim.complex_fft is True else ref + sw / 2.0

        if sw < 0:
            sw = -sw
            co = ref + increment * ((count + 1) // 2 - 1)

        dim_i = {}
        dim_i["count"] = dim.count
        dim_i["spectral_width"] = sw
        dim_i["reference_offset"] = co

        if dim.label != "":
            dim_i["label"] = dim.label
        if dim.origin_offset not in ["", None, 0]:
            dim_i["origin_offset"] = dim.origin_offset.to("Hz").value
        result.append(dim_i)

    return result[::-1]


def flatten_dict(obj, previous_key=None):
    """Flatten a nested dictionary with keys as obj1.obj2... and so on"""
    result = {}
    for k, v in obj.items():
        if not isinstance(v, dict):
            key = f"{previous_key}.{k}" if previous_key is not None else k
            result.update({key: v})
        else:
            result.update(**flatten_dict(v, previous_key=k))

    return result
