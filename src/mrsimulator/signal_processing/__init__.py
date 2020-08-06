# -*- coding: utf-8 -*-
"""The Event class."""
from sys import modules
from typing import List
from typing import Union

import csdmpy as cp
from pydantic import BaseModel

from ._base import AbstractOperation

__author__ = "Maxwell C. Venetos"
__email__ = "maxvenetos@gmail.com"


class SignalProcessor(BaseModel):
    """
    Signal processing class to apply a series of operations to the dependent variables
    of the simulation dataset.

    Attributes
    ----------

    operations: List
        A list of operations.

    Examples
    --------

    >>> post_sim = SignalProcessor(operations=[o1, o2]) # doctest: +SKIP
    """

    processed_data: cp.CSDM = None
    operations: List[AbstractOperation] = []

    class Config:
        validate_assignment = True
        arbitrary_types_allowed = True

    @classmethod
    def parse_dict_with_units(self, py_dict):
        """Parse a list of operations dictionary to a SignalProcessor class object.

        Args:
            pt_dict: A python dict object.
        """
        lst = []
        for op in py_dict["operations"]:
            if op["function"] == "apodization":
                lst.append(
                    getattr(
                        modules[__name__].apodization, op["type"]
                    ).parse_dict_with_units(op)
                )
            else:
                lst.append(
                    getattr(modules[__name__], op["function"]).parse_dict_with_units(op)
                )
        return SignalProcessor(operations=lst)

    def to_dict_with_units(self):
        """
        Serialize the SignalProcessor object to a JSON compliant python dictionary
        object, where physical quantities are represented as string with a value and a
        unit.

        Returns:
            A Dict object.
        """
        lst = []
        for i in self.operations:
            lst += [i.to_dict_with_units()]
        op = {}

        op["operations"] = lst
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
        copy_data = data.copy()
        for filters in self.operations:
            copy_data = filters.operate(copy_data)
        self.processed_data = copy_data

        return copy_data


class Scale(AbstractOperation):
    """
    Scale the amplitudes of all dependent variables from a CSDM object.

    Args:
        float factor: The scaling factor. The default value is 1.

    Example
    -------

    >>> import mrsimulator.signal_processing as sp
    >>> operation1 = sp.Scale(factor=20)
    """

    factor: float = 1

    def operate(self, data):
        r"""Applies the operation for which the class is named for.

        .. math::
            f(\vec(x)) = scale*\vec(x)

        Args:
            data: CSDM object
        """
        data *= self.factor
        return data


class IFFT(AbstractOperation):
    """
    Apply an inverse Fourier transform on all dependent variables of the CSDM object.

    Args:
        int dim_index: Dimension index along which the function is applied.

    Example
    -------

    >>> operation2 = sp.IFFT(dim_index=0)
    """

    dim_index: Union[int, list, tuple] = 0

    def operate(self, data):
        """Applies the operation for which the class is named for.

        Args:
            data: CSDM object
        """
        dim_index = self.dim_index
        if isinstance(dim_index, int):
            dim_index = [dim_index]

        for i in dim_index:
            data = data.fft(axis=i)
        return data


class FFT(IFFT):
    """
    Apply a forward Fourier transform on all dependent variables of the CSDM object.

    Args:
        int dim_index: Dimension index along which the function is applied.

    Example
    -------

    >>> operation3 = sp.FFT(dim_index=0)
    """

    pass


class complex_conjugate(AbstractOperation):
    pass
