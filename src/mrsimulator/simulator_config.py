# -*- coding: utf-8 -*-
"""Base ConfigSimulator class."""
from copy import deepcopy

from mrsimulator.sandbox import AveragingScheme


__author__ = "Deepansh J. Srivastava"
__email__ = "deepansh2012@gmail.com"


__integration_volume_enum__ = {"octant": 0, "hemisphere": 1}
__integration_volume_enum_rev__ = {0: "octant", 1: "hemisphere"}


class ConfigSimulator:
    """
    The configurable parametes used in lineshape simulation

    Attributes:
        number_of_sidebands: An integer with the number of sidebands requested in the
            lineshape simulation. This number cannot be zero or negative. The default
            value is ``64``.
        integration_volume: An enumeration literal with valid literals, 'octant',
            'hemisphere', over which the frequency contributions are evaluated. The
            default value if ``octant``.
        integration_density: An positive interger. If `n` is the integration_density,
            then the total number of orientation is given as
            ``(n+1)*(n+2)/2 * number of octant``. The line-shape is an integral over
            the frequency contribution arising from every orientation. The default
            value is ``70``.
        decompose: A boolean. If true, decomposes the line-shape into an array of
            line-shapes arising from individual isotopomer. If False, the lins-shape
            is a sum of individual line-shapes instead. The default value is ``False``.

    Example:
        >>> a = Simulator()
        >>> a.config.number_of_sidebands = 128
        >>> a.config.integration_volume = 'hemisphere'
        >>> a.config.decompose = True
    """

    def __init__(self):
        self._dict = {
            "number_of_sidebands": 64,
            "integration_volume": 0,
            "integration_density": 70,
            "decompose": False,
        }
        self._averaging_scheme = AveragingScheme(
            integration_density=70, integration_volume=0, allow_fourth_rank=True
        )

    # def update_averaging_scheme(self):
    #     self._averaging_scheme = AveragingScheme(
    #         integration_density=self._dict["integration_density"],
    #         integration_volume=self._dict["integration_volume"],
    #         allow_fourth_rank=True,
    #     )

    # decompose line-shape into array of lineshapes arising from individual
    # isotopomers.
    @property
    def decompose(self):
        return self._dict["decompose"]

    @decompose.setter
    def decompose(self, value):
        if isinstance(value, bool):
            self._dict["decompose"] = value
            return
        raise ValueError(f"Expecting a boolean.")

    # number of sidebands
    @property
    def number_of_sidebands(self):
        return self._dict["number_of_sidebands"]

    @number_of_sidebands.setter
    def number_of_sidebands(self, value):
        if isinstance(value, int):
            if value >= 0:
                self._dict["number_of_sidebands"] = value
                return
        raise ValueError(f"Expecting a positive integer, found {value}.")

    # integration density
    @property
    def integration_density(self):
        # return self._averaging_scheme.integration_density
        return self._dict["integration_density"]

    @integration_density.setter
    def integration_density(self, value):
        # self._averaging_scheme.integration_density = value
        # self._dict["integration_density"] = value
        if isinstance(value, int):
            if value > 0:
                self._dict["integration_density"] = value
                # self.update_integration_scheme()
                return
        raise ValueError(f"Expecting a positive integer, found {value}.")

    # integration volume
    @property
    def integration_volume(self):
        # return self._averaging_scheme.integration_volume
        return __integration_volume_enum_rev__[self._dict["integration_volume"]]

    @integration_volume.setter
    def integration_volume(self, value):
        # self._averaging_scheme.integration_volume = value
        # self._dict["integration_volume"] = value
        # class_dict = self.__class__.integration_volume_enum
        if value in __integration_volume_enum__.keys():
            self._dict["integration_volume"] = __integration_volume_enum__[value]
            return
        raise ValueError(
            (
                "value is not a valid enumeration literal; "
                "permitted: 'octant', 'hemisphere', found {value}.",
            )
        )

    # averaging scheme. This contains the c pointer used in line-shape evaluation
    # at c-level.
    @property
    def averaging_scheme(self):
        return self._averaging_scheme

    @averaging_scheme.setter
    def averaging_scheme(self, other):
        if isinstance(other, AveragingScheme):
            self._averaging_scheme = other
        raise ValueError("Expecting an instance of either the AveragingScheme class.")

    # methods
    def __eq__(self, other):
        if isinstance(other, ConfigSimulator):
            if self._dict == other._dict:
                return True
        return False

    def __str__(self):
        return str(self.dict())

    def __repr__(self):
        return (
            "ConfigSimulator(number_of_sidebands={0}, integration_volume={1}, "
            "integration_density={2}, decompose={3})"
        ).format(
            # self.averaging_scheme,
            self.number_of_sidebands,
            self.integration_volume,
            self.integration_density,
            self.decompose,
        )

    # dict
    def dict(self):
        dict_self = deepcopy(self._dict)
        dict_self["integration_volume"] = self.integration_volume
        return dict_self
