# -*- coding: utf-8 -*-
from sys import modules
from typing import List
from typing import Union

import csdmpy as cp
from mrsimulator.utils.parseable import Parseable
from pydantic import Extra

from . import affine as af  # noqa:F401
from . import apodization as ap  # noqa:F401
from ._base import Operations

__author__ = "Maxwell C. Venetos"
__email__ = "maxvenetos@gmail.com"


class SignalProcessor(Parseable):
    """
    Signal processing class to apply a series of operations to the dependent variables
    of the simulation dataset.

    Arguments
    ---------

    operations:
        A list of operations.

    Examples
    --------

    >>> post_sim = SignalProcessor(operations=[o1, o2]) # doctest: +SKIP
    """

    processed_data: cp.CSDM = None
    operations: List[Operations] = []

    class Config:
        validate_assignment = True
        arbitrary_types_allowed = True
        extra = Extra.forbid

    @classmethod
    def parse_dict_with_units(cls, py_dict: dict):
        """Parse a list of operations dictionary to a SignalProcessor class object.

        Args:
            pt_dict: A python dict object.
        """
        mod = modules[__name__]
        lst = [
            getattr(getattr(mod, op["function"]), op["type"]).parse_dict_with_units(op)
            if "type" in op
            else getattr(mod, op["function"]).parse_dict_with_units(op)
            for op in py_dict["operations"]
        ]
        return SignalProcessor(operations=lst)

    def json(self) -> dict:
        """Parse the class object to a JSON compliant python dictionary object, where
        the attribute value with physical quantity is expressed as a string with a
        value and a unit.

        Returns:
            A Dict object.
        """
        op = {}
        op["operations"] = [item.json() for item in self.operations]
        return op

    def apply_operations(self, data, **kwargs):
        """
        Function to apply all the operation functions in the operations member of a
        SignalProcessor object. Operations applied sequentially over the data member.

        Returns:
            CSDM object: A copy of the data member with the operations applied to it.
        """
        if not isinstance(data, cp.CSDM):
            raise ValueError("The data must be a CSDM object.")
        for filters in self.operations:
            data = filters.operate(data)
        self.processed_data = data

        return data


class Scale(Operations):
    r"""
    Scale the amplitudes of all dependent variables from a CSDM object.

    .. math::
        f(\vec(x)) = scale*\vec(x)

    Arguments
    ---------

    factor:
        The scaling factor. The default value is 1.

    Example
    -------

    >>> from mrsimulator import signal_processing as sp
    >>> operation1 = sp.Scale(factor=20)
    """

    factor: float = 1
    function: str = "Scale"

    def operate(self, data):
        r"""Applies the operation.
        Args:
            data: CSDM object
        """
        data *= self.factor
        return data


class IFFT(Operations):
    """
    Apply an inverse Fourier transform on all dependent variables of the CSDM object.

    Arguments
    ---------

    dim_index:
        Dimension index along which the function is applied.

    Example
    -------

    >>> operation2 = sp.IFFT(dim_index=0)
    """

    dim_index: Union[int, list, tuple] = 0
    function: str = "IFFT"

    def operate(self, data):
        r"""Applies the operation.
        Args:
            data: CSDM object
        """
        d_i = [self.dim_index] if isinstance(self.dim_index, int) else self.dim_index
        for i in d_i:
            data = data.fft(axis=i)
        return data


class FFT(IFFT):
    """
    Apply a forward Fourier transform on all dependent variables of the CSDM object.

    Arguments
    ---------

    dim_index:
        Dimension index along which the function is applied.

    Example
    -------

    >>> operation3 = sp.FFT(dim_index=0)
    """

    function: str = "FFT"


class complex_conjugate(Operations):
    pass
