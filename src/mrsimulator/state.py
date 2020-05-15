# -*- coding: utf-8 -*-
class ZeemanState:
    """Zeeman energy state class."""

    def __init__(self, n_sites, *args):
        self.n_sites = n_sites
        for i in range(n_sites):
            self.__setattr__(f"m{i}", args[i])

    def __repr__(self):
        return self.__str__()
        # lst = "".join(
        #     [f"m{i}={self.__getattribute__(f'm{i}')}, " for i in range(self.n_sites)]
        # )
        # return f"ZeemanState({lst[:-2]})"

    def __str__(self):
        lst = "".join(
            [f"{self.__getattribute__(f'm{i}')}, " for i in range(self.n_sites)]
        )
        return f"|{lst[:-2]}⟩"

        # lst = "".join(
        #     [f"m_{i}={self.__getattribute__(f'm_{i}')} " for i in range(self.n_sites)]
        # )
        # return f"{lst[:-1]}"

    def tolist(self):
        return [self.__getattribute__(f"m{i}") for i in range(self.n_sites)]
