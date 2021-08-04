# -*- coding: utf-8 -*-

__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"


class ZeemanState:
    """Zeeman energy state class."""

    def __init__(self, n_sites, *args):
        self.n_sites = n_sites
        for i in range(n_sites):
            self.__setattr__(f"m{i}", args[i])

    def __repr__(self):
        return self.__str__()

    # def _repr_html_(self):
    #   lst = [int(2 * self.__getattribute__(f"m{i}")) for i in range(self.n_sites)]
    #   string = " ".join([rf"m_{i}={{{item} \over 2}}" for i, item in enumerate(lst)])
    #   return rf"$|{string}⟩$"

    def __str__(self):
        lst = ", ".join(
            [f"{self.__getattribute__(f'm{i}')}" for i in range(self.n_sites)]
        )
        return f"|{lst}⟩"

    def tolist(self):
        return [self.__getattribute__(f"m{i}") for i in range(self.n_sites)]
