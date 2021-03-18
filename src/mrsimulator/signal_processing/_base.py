# -*- coding: utf-8 -*-
"""The Operations class."""
import numpy as np
from mrsimulator.utils.parseable import Parseable

__author__ = "Maxwell C. Venetos"
__email__ = "maxvenetos@gmail.com"


class Operations(Parseable):
    """A base class for signal processing operations."""

    @staticmethod
    def _get_dv_indexes(indexes, n):
        """Return a list of dependent variable indexes.

        Args:
            indexes: An interger, list of integers, or None indicating the dv indexes.
            n: Total number of dependent variables in the CSDM object.
        """
        if indexes is None:
            return np.arange(n)
        if isinstance(indexes, int):
            return [indexes]
        if isinstance(indexes, (list, tuple)):
            return np.asarray(indexes)
