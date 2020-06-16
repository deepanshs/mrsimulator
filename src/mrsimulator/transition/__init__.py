# -*- coding: utf-8 -*-
"""The Transition class."""
from typing import List

import numpy as np
from pydantic import BaseModel

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


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

    def __repr__(self):
        """Representation in bar-ket notation, |final⟩⟨initial|."""
        final = ", ".join([str(i) for i in self.final])
        initial = ", ".join([str(i) for i in self.initial])
        return f"|{final}⟩⟨{initial}|"

    # @property
    # def Zeeman_allowed(self):
    #     if abs(self.p) == 1:
    #         return True
    #     return False

    @property
    def p(self):
        """Return the total Δm (m_final-m_initial) value of the spin transition."""
        return self.P.sum()

    @property
    def delta_m(self):
        """An alias for p"""
        return self.p

    @property
    def P(self):
        """Return a list of Δm values of the spin transition for each site."""
        return np.asarray(self.final) - np.asarray(self.initial)

    @property
    def D(self):
        """Return a list of Δm**2 values of the spin transition for each site."""
        return np.asarray(self.final) ** 2 - np.asarray(self.initial) ** 2

    # def delta_m(self, i):
    #     """Return the Δm element of the transition corresponding to the ith site."""
    #     return self.final[i] - self.initial[i]

    def tolist(self):
        """Convert the transition to a list of quantum numbers where the first N
        quantum numbers corresponds to the initial energy state while the last N
        corresponds to the final energy state, where N is the number of sites."""
        lst = self.initial + self.final
        return lst

    def to_dict_with_units(self):
        return self.dict()
