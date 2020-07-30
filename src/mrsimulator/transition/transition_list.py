# -*- coding: utf-8 -*-
import numpy as np
from mrsimulator.utils.abstract_list import AbstractList

from . import Transition

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


class TransitionList(AbstractList):
    def __init__(self, data=[]):
        super().__init__([self._check_for_transition_object(item) for item in data])

    @staticmethod
    def _check_for_transition_object(item):
        if isinstance(item, dict):
            return Transition(**item)
        if not isinstance(item, Transition):
            raise ValueError(
                "Expecting a Transition object or an equivalent python dict object, "
                f"instead found {item.__class__.__name__}."
            )
        return item

    def __setitem__(self, index, item):
        self._list[index] = self._check_for_transition_object(item)

    def append(self, item):
        super().append(self._check_for_transition_object(item))

    # def Zeeman_allowed(self):
    #     return TransitionList([item for item in self._list if item.Zeeman_allowed])

    def filter(self, P=None, D=None):
        """
        Filter a list of transitions to satisfy the filtering criterion.
        Args:
            list P: A list of `N` (m_final - m_initial) values, where `N` is the
                total number of sites within the spin system.
            list D: A list of `N` (m_final^2 - m_initial^2) values, where `N` is the
                total number of sites within the spin system.
        """

        # to think
        # - filter based on transition.
        # - filter based on state.

        if P is D is None:
            return self

        ts = self._list.copy()

        if P is not None:
            ts = TransitionList([item for item in ts if np.allclose(item.P, P)])
        if D is not None:
            ts = TransitionList([item for item in ts if np.allclose(item.D, D)])
        # if transitions is not None:
        #     for transition in transitions:
        #         ts = TransitionList(
        #             [item for item in ts if item == Transition(**transition)]
        #         )
        # if start_state is not None:
        #     ts = [item for item in ts if ts.initial == start_state]
        return ts

        # lst = self._list
        # for i, element in enumerate(search):
        #     if isinstance(element, int):
        #         lst = [item for item in lst if item.delta_m(i) == element]
        #     if isinstance(element, list):
        #         lst = [
        #             item
        #             for item in lst
        #             if (item.initial[i] == element[1] and item.final[i] == element[0])
        #         ]
        # return lst
        #     return TransitionList(
        #         [item for item in self._list if np.all(item.initial == delta_ms)]
        #     )

        # if delta_ms is not None:
        #     return TransitionList(
        #         [item for item in self._list if np.all(item.delta_ms == delta_ms)]
        #     )
        # return TransitionList([item for item in self._list if item.p == p])
