# -*- coding: utf-8 -*-
"""The Transition class."""
from typing import List

import numpy as np
from pydantic import BaseModel


class Transition(BaseModel):
    """
    Base Transition class.

    Arguments:
        initial: The initial Zeeman energy state represented as a list of quantum
            numbers.
        final: The final Zeeman energy state represented as a list of quantum
            numbers.
    """

    initial: List[float] = []
    final: List[float] = []

    @property
    def Zeeman_allowed(self):
        if abs(self.p) == 1:
            return True
        return False

    @property
    def p(self):
        """Return the total Δm value of the spin transition."""
        return self.delta_ms.sum()

    @property
    def delta_ms(self):
        """Return a list of Δm values of the spin transition for each site."""
        return np.asarray(self.final) - np.asarray(self.initial)

    def delta_m(self, i):
        """Return the Δm element of the transition corresponding to the ith site."""
        return self.final[i] - self.initial[i]
