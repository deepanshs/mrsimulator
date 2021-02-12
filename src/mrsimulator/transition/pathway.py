# -*- coding: utf-8 -*-
import numpy as np
from mrsimulator.utils.abstract_list import AbstractList

from . import Transition

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


class TransitionList(AbstractList):
    def __init__(self, data: list = []):
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

    def filter(self, P=None, PP=None, D=None):
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

        if P is PP is D is None:
            return self

        ts = self._list.copy()

        if P is not None:
            ts = TransitionList([item for item in ts if np.allclose(item.P, P)])

        # if PP is not None:
        #     ts = TransitionList([item for item in ts if np.allclose(item.PP, PP)])

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


class TransitionPathway(TransitionList):
    """
    Base TransitionPathway class is a list of connected Transitions.

    Example:
        >>> from mrsimulator.transition import TransitionPathway, Transition
        >>> t1 = Transition(initial = [0.5, 0.5], final = [0.5, -0.5])
        >>> t2 = Transition(initial=[0.5, 0.5], final=[-0.5, 0.5])
        >>> path = TransitionPathway([t1, t2])
        >>> path
        |0.5, -0.5⟩⟨0.5, 0.5| ⟶ |-0.5, 0.5⟩⟨0.5, 0.5|
    """

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return " ⟶ ".join([repr(item) for item in self._list])

    def json(self) -> dict:
        """Parse the class object to a JSON compliant python dictionary object.

        Example:
            >>> pprint(path.json())
            [{'final': [0.5, -0.5], 'initial': [0.5, 0.5]},
             {'final': [-0.5, 0.5], 'initial': [0.5, 0.5]}]
        """
        return [item.json() for item in self._list]

    def tolist(self):
        """Expand TransitionPathway to a Python list.

        Example:
            >>> path.tolist()
            [0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5]
        """
        return list(np.asarray([item.tolist() for item in self._list]).ravel())
