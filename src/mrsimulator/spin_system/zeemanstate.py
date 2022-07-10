__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"


class ZeemanState:
    """Zeeman energy state class."""

    def __init__(self, n_sites, *args):
        self.n_sites = n_sites
        _ = [setattr(self, f"m{i}", args[i]) for i in range(n_sites)]

    def __repr__(self):
        return self.__str__()

    # def _repr_html_(self):
    #   lst = [int(2 * getattr(self, f"m{i}")) for i in range(self.n_sites)]
    #   string = " ".join([rf"m_{i}={{{item} \over 2}}" for i, item in enumerate(lst)])
    #   return rf"$|{string}⟩$"

    def __str__(self):
        lst = ", ".join([f"{getattr(self, f'm{i}')}" for i in range(self.n_sites)])
        return f"|{lst}⟩"

    def tolist(self):
        return [getattr(self, f"m{i}") for i in range(self.n_sites)]
