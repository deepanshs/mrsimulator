# -*- coding: utf-8 -*-
"""The Transition class."""
from typing import List

import numpy as np
from pydantic import BaseModel

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


class Transition(BaseModel):
    r"""
    Base Transition class describes a spin transition between two energy states, where
    the energy states are described using the weakly coupled basis.

    .. math::
        |m_{i,0}, m_{i,1}, ... m_{i,N} \rangle \rightarrow
            |m_{f,0}, m_{f,1}, ... m_{f,N} \rangle

    Arguments:
        list initial: The initial Zeeman energy state represented as a list of quantum
            numbers :math:`m_{i,n}`.
        list final: The final Zeeman energy state represented as a list of quantum
            numbers :math:`m_{f,n}`.

    Example:
        >>> from mrsimulator.transition import Transition
        >>> t1 = Transition(initial = [0.5, 0.5], final = [0.5, -0.5])
        >>> t1
        |0.5, -0.5⟩⟨0.5, 0.5|
    """

    initial: List[float] = []
    final: List[float] = []

    def __repr__(self):
        """Representation in bar-ket notation, |final⟩⟨initial|."""
        final = ", ".join([str(i) for i in self.final])
        initial = ", ".join([str(i) for i in self.initial])
        return f"|{final}⟩⟨{initial}|"

    def __str__(self):
        return self.__repr__()

    # @property
    # def Zeeman_allowed(self):
    #     if abs(self.p) == 1:
    #         return True
    #     return False

    @property
    def p(self):
        """Return the total Δm (m_final-m_initial) value of the spin transition.

        Example:
            >>> t1.p
            -1.0
        """
        return self.P.sum()

    @property
    def delta_m(self):
        """An alias for p"""
        return self.p

    @property
    def P(self):
        """Return a list of Δm values of the spin transition for each site in a
        weakly coupled basis.

        Example:
            >>> t1.P
            array([ 0., -1.])
        """
        return np.asarray(self.final) - np.asarray(self.initial)

    # @property
    # def PP(self):
    #     """Return a list of Δm values of the spin transition for each site."""
    #     return np.prod(np.asarray(self.final)) - np.prod(np.asarray(self.initial))

    @property
    def D(self):
        """Return a list of Δm**2 values of the spin transition for each site in a
        weakly coupled basis.

        Example:
            >>> t1.D
            array([0., 0.])
        """
        return np.asarray(self.final) ** 2 - np.asarray(self.initial) ** 2

    # def delta_m(self, i):
    #     """Return the Δm element of the transition corresponding to the ith site."""
    #     return self.final[i] - self.initial[i]

    def tolist(self) -> list:
        """Convert the transition to a list of quantum numbers where the first N
        quantum numbers corresponds to the initial energy state, while the last N
        corresponds to the final energy state, where N is the number of sites.

        Example:
            >>> t1.tolist()
            [0.5, 0.5, 0.5, -0.5]
        """
        lst = self.initial + self.final
        return lst

    def json(self) -> dict:
        """Parse the class object to a JSON compliant python dictionary object.

        Example:
            >>> t1.json()
            {'initial': [0.5, 0.5], 'final': [0.5, -0.5]}
        """
        return self.dict()
