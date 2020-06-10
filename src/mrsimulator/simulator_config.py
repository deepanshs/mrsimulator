# -*- coding: utf-8 -*-
"""Base ConfigSimulator class."""
from copy import deepcopy

from mrsimulator.sandbox import AveragingScheme


__author__ = "Deepansh J. Srivastava"
__email__ = "deepansh2012@gmail.com"

# decompose spectrum
__decompose_spectrum_enum__ = {"none": 0, "spin_system": 1}
__decompose_spectrum_enum_rev__ = {0: "none", 1: "spin_system"}

# integration volume
__integration_volume_enum__ = {"octant": 0, "hemisphere": 1}
__integration_volume_enum_rev__ = {0: "octant", 1: "hemisphere"}
__integration_volume_octants__ = [1, 4]


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
        decompose_spectrum: A boolean. If true, decomposes the line-shape into an array
            of line-shapes arising from individual isotopomer. If False, the lins-shape
            is a sum of individual line-shapes instead. The default value is ``False``.

    Example
    -------

    >>> a = Simulator()
    >>> a.config.number_of_sidebands = 128
    >>> a.config.integration_volume = 'hemisphere'
    >>> a.config.decompose_spectrum = 'spin_system'
    """

    def __init__(self, **kwargs):
        self._dict = {
            "number_of_sidebands": 64,
            "integration_volume": 0,
            "integration_density": 70,
            "decompose_spectrum": 0,
        }
        if kwargs != {}:
            self._dict = {
                "number_of_sidebands": kwargs["number_of_sidebands"],
                "integration_volume": __integration_volume_enum__[
                    kwargs["integration_volume"]
                ],
                "integration_density": kwargs["integration_density"],
                "decompose_spectrum": __decompose_spectrum_enum__[
                    kwargs["decompose_spectrum"]
                ],
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
    # spin systems.
    @property
    def decompose_spectrum(self):
        return __decompose_spectrum_enum_rev__[self._dict["decompose_spectrum"]]

    @decompose_spectrum.setter
    def decompose_spectrum(self, value):
        if not isinstance(value, str):
            raise TypeError("Expecting a string value.")
        if value in __decompose_spectrum_enum__.keys():
            self._dict["decompose_spectrum"] = __decompose_spectrum_enum__[value]
            return
        raise ValueError(
            (
                "Value is not a valid enumeration literal; "
                "permitted: 'none', 'spin_system', found {value}.",
            )
        )

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
        if not isinstance(value, str):
            raise TypeError("Expecting a string value.")
        if value in __integration_volume_enum__.keys():
            self._dict["integration_volume"] = __integration_volume_enum__[value]
            return
        raise ValueError(
            (
                "Value is not a valid enumeration literal; "
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
            "integration_density={2}, decompose_spectrum={3})"
        ).format(
            # self.averaging_scheme,
            self.number_of_sidebands,
            self.integration_volume,
            self.integration_density,
            self.decompose_spectrum,
        )

    # dict
    def dict(self):
        dict_self = deepcopy(self._dict)
        dict_self["integration_volume"] = self.integration_volume
        dict_self["decompose_spectrum"] = self.decompose_spectrum
        return dict_self

    def get_orientations_count(self):
        """Return the total number of orientations."""
        n = self._dict["integration_density"]
        vol = __integration_volume_octants__[self._dict["integration_volume"]]
        return int(vol * (n + 1) * (n + 2) / 2)
