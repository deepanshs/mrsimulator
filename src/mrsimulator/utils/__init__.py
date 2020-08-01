# -*- coding: utf-8 -*-


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
        coordinates_offset = dim.coordinates_offset.to("Hz").value
        spectral_width = count * increment

        dim_i = {}
        dim_i["count"] = dim.count
        dim_i["spectral_width"] = spectral_width
        if dim.complex_fft is True:
            dim_i["reference_offset"] = coordinates_offset
        else:
            dim_i["reference_offset"] = coordinates_offset + spectral_width / 2.0

        result.append(dim_i)

    return result


def _reduce_dict(dict_obj, exclude=["property_units"]):
    """Reduce the dict by removing all key-value pair corresponding to keys listed in
    the `exclude` argument and keys with value as None.

    Args:
        exclude: List of keys to exclude from the dict.
    Returns: A dict.
    """
    if isinstance(dict_obj, dict):
        obj = {}
        for k, v in dict_obj.items():
            if k in exclude or v is None:
                continue
            if k == "isotope":
                obj[k] = v["symbol"]
                continue
            if k == "channels":
                obj[k] = [item["symbol"] for item in v]
                continue
            elif isinstance(v, dict):
                obj[k] = _reduce_dict(v)
            elif isinstance(v, list):
                obj[k] = [_reduce_dict(_) for _ in v]
            else:
                obj[k] = v
        return obj

    # if isinstance(dict_obj, list):
    #     obj = []
    #     for v in dict_obj:
    #         if v is None:
    #             continue
    #         if isinstance(v, dict):
    #             obj.append(_reduce_dict(v))
    #         elif isinstance(v, list):
    #             obj.append([_reduce_dict(_) for _ in v])
    #         else:
    #             obj.append(v)
    #     return obj

    return dict_obj
