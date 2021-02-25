# -*- coding: utf-8 -*-
from mrsimulator import methods as NamedMethods


__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"

methods = [val for k, val in NamedMethods.__dict__.items() if isinstance(val, type)]

sp = {"count": 1024, "spectral_width": 100}


def test_read_write_methods():
    for item in methods:
        kwargs = {"channels": ["93Nb"], "spectral_dimensions": [sp] * item.ndim}

        if item.__name__ == "SSB2D":
            kwargs["rotor_frequency"] = 1e3

        fn1 = item(**kwargs)
        fn2 = item.parse_dict_with_units(fn1.json())
        assert fn1 == fn2, f"Error with {item} parse."
