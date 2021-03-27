# -*- coding: utf-8 -*-

IGNORE = ["simulation", "experiment"]


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
            if k in IGNORE and v is not None:
                obj[k] = v
                continue
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
