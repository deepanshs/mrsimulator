# -*- coding: utf-8 -*-
"""The Event class."""
from typing import List

import csdmpy as cp
from pydantic import BaseModel

from ._base import abstractOperation

__author__ = "Maxwell C. Venetos"
__email__ = "maxvenetos@gmail.com"


class SignalProcessor(BaseModel):
    """
    Signal processing class to apply lists of various operations to individual dependent variables
    of the data.

    Attributes
    ----------

    data: CSDM object.
        From simulation

    operations: List
        List of operation lists

    Examples
    --------

    >>> post_sim = SignalProcessor(data = csdm_object, operations = [op1, op2]) # doctest: +SKIP
    """

    data: cp.CSDM = None
    operations: List[abstractOperation]

    class Config:
        validate_assignment = True
        arbitrary_types_allowed = True

    def to_dict_with_units(self):
        """
        Serialize the SignalProcessor object to a JSON compliant python dictionary object
        where physical quantities are represented as string with a value and a unit.

        Returns:
            A Dict object.
        """
        lst = []
        for i in self.operations:
            lst += [i.to_dict_with_units()]
        op = {}

        op["operations"] = lst
        return op

    def apply_operations(self, **kwargs):
        """
        Function to apply all the operation functions in the operations member of a
        SignalProcessor object. Operations applied sequentially over the data member.

        Returns:
            CSDM object: A copy of the data member with the operations applied to it.
        """
        copy_data = self.data.copy()
        for filters in self.operations:
            copy_data = filters.operate(copy_data)

        return copy_data


class Scale(abstractOperation):
    """
    Class for applying a scaling
    factor to a dependent variable of simulation data
    """

    factor: float = 1

    def operate(self, data):
        r"""
        Applies the operation for which the class is named for.

        .. math::
            f(\vec(x)) = scale*\vec(x)

        Args:
            data: CSDM object
        """
        data *= self.factor
        return data


class IFFT(abstractOperation):
    """
    Class for applying an inverse
    Fourier transform to a dependent variable of simulation data

    Args:
        dim_indx: int. Data dimension to apply the function along

    """

    dim_indx: int = 0

    def operate(self, data):
        """
        Applies the operation for which the class is named for.

        Args:
            data: CSDM object

        """
        return data.fft(axis=self.dim_indx)


class FFT(IFFT):
    pass


class complex_conjugate(abstractOperation):
    pass
