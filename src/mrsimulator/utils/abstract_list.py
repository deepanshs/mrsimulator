# -*- coding: utf-8 -*-
from collections.abc import MutableSequence

import numpy as np

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


class AbstractList(MutableSequence):
    def __init__(self, data=[]):
        super().__init__()
        self._list = list(data)

    def __repr__(self):
        """String representation"""
        return self.__str__()

    def __str__(self):
        """String representation"""
        string = ",\n".join([item.__repr__() for item in self._list])
        return f"[{string}]"

    def __len__(self):
        """List length"""
        return len(self._list)

    def __getitem__(self, index):
        """Get a list item"""
        return self._list[index]

    def __delitem__(self, index):
        """Delete an item from the list."""
        del self._list[index]

    def insert(self, index, item):
        """Insert a list item"""
        self._list.insert(index, item)

    def append(self, item):
        """Append a list item"""
        self.insert(len(self._list), item)

    def __eq__(self, other):
        """Check equality of DependentVariableList."""
        if not isinstance(other, self.__class__):
            return False
        if len(self._list) != len(other._list):
            return False

        check = []
        for i in range(len(self._list)):
            check.append(self._list[i] == other._list[i])

        if np.all(check):
            return True
        return False
