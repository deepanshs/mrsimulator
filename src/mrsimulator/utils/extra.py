import numpy as np


__author__ = ["Matthew D. Giammar"]
__email__ = ["giammar.7@buckeyemail.osu.edu"]


def broadcast_dict_values(d, n_sites=None):
    """Checks the type of values in the passed dictionary, d, and extends singular
    values 1d numpy arrays of length n_sites.

    Example:
    >>> d = {"foo": 1, "bar": [6, 7, 8, 9]}
    >>> broadcast_dict_values(d, 4)
    {'foo': array([1, 1, 1, 1]), 'bar': array([6, 7, 8, 9])}
    """
    n_sites_dict = _check_dict_value_lengths(d)
    if n_sites is None:
        n_sites = n_sites_dict
    elif n_sites_dict != n_sites:
        raise ValueError("Length of lists in passed dict did not match `n_sites`")

    return {
        key: np.ravel(np.asarray(val))
        if isinstance(val, (list, np.ndarray))
        else np.full((n_sites), val)
        for key, val in d.items()
    }


def zip_dict(d):
    """Returns a list of dicts when the dict passed has values which are lists or
    numpy arrays. Singular values will be broadcast to all dictionaries

    Example:
    >>> d = {"foo": 1, "bar": [7, 8, 9]}
    >>> zip_dict(d)
    [{'foo': 1, 'bar': 7}, {'foo': 1, 'bar': 8}, {'foo': 1, 'bar': 9}]
    """
    d = broadcast_dict_values(d)
    lst = [dict(zip(d.keys(), v)) for v in zip(*(d[k] for k in d.keys()))]
    return lst


def _check_dict_value_lengths(d):
    """Raises value error when a dict has values of multiple lengths"""
    length = {len(v) if isinstance(v, (list, np.ndarray)) else None for v in d.values()}
    length.discard(None)

    if len(length) == 0:
        return 1

    # If different lengths present in dictionary, throw error
    if len(length) > 1:
        raise ValueError(
            f"Values with multiple lengths passed, {length}. "
            "Can't infer a length to broadcast with values of multiple lengths."
        )
    return length.pop()


# def _reduce_dict(dict_obj, exclude=["property_units"]):
#     """Reduce the dict by removing all key-value pair corresponding to keys listed in
#     the `exclude` argument and keys with value as None.
#     Args:
#         exclude: List of keys to exclude from the dict.
#     Returns: A dict.
#     """
#     if isinstance(dict_obj, dict):
#         obj = {}
#         for k, v in dict_obj.items():
#             if k in exclude or v is None:
#                 continue
#             if k == "isotope":
#                 obj[k] = v["symbol"]
#                 continue
#             if k == "channels":
#                 obj[k] = [item["symbol"] for item in v]
#                 continue
#             elif isinstance(v, dict):
#                 obj[k] = _reduce_dict(v)
#             elif isinstance(v, list):
#                 obj[k] = [_reduce_dict(_) for _ in v]
#             else:
#                 obj[k] = v
#         return obj
#     # if isinstance(dict_obj, list):
#     #     obj = []
#     #     for v in dict_obj:
#     #         if v is None:
#     #             continue
#     #         if isinstance(v, dict):
#     #             obj.append(_reduce_dict(v))
#     #         elif isinstance(v, list):
#     #             obj.append([_reduce_dict(_) for _ in v])
#     #         else:
#     #             obj.append(v)
#     #     return obj
#     return dict_obj
