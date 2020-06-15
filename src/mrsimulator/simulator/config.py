# -*- coding: utf-8 -*-
"""Base ConfigSimulator class."""
from mrsimulator.sandbox import AveragingScheme
from pydantic import BaseModel
from pydantic import Field
from typing_extensions import Literal

__author__ = "Deepansh J. Srivastava"
__email__ = "deepansh2012@gmail.com"

# decompose spectrum
__decompose_spectrum_enum__ = {"none": 0, "spin_system": 1}

# integration volume
__integration_volume_enum__ = {"octant": 0, "hemisphere": 1}
__integration_volume_octants__ = [1, 4]


class ConfigSimulator(BaseModel):
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
            of line-shapes arising from individual spin systems. If False, the lins-shape
            is a sum of individual line-shapes instead. The default value is ``False``.

    Example
    -------

    >>> a = Simulator()
    >>> a.config.number_of_sidebands = 128
    >>> a.config.integration_volume = 'hemisphere'
    >>> a.config.decompose_spectrum = 'spin_system'
    """

    number_of_sidebands: int = Field(default=64, gt=0)
    integration_volume: Literal["octant", "hemisphere"] = "octant"
    integration_density: int = Field(default=70, gt=0)
    decompose_spectrum: Literal["none", "spin_system"] = "none"

    class Config:
        validate_assignment = True

    def get_int_dict(self):
        py_dict = self.dict()
        py_dict["integration_volume"] = __integration_volume_enum__[
            self.integration_volume
        ]
        py_dict["decompose_spectrum"] = __decompose_spectrum_enum__[
            self.decompose_spectrum
        ]
        return py_dict

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

    def get_orientations_count(self):
        """Return the total number of orientations."""
        n = self.integration_density
        vol = __integration_volume_octants__[
            __integration_volume_enum__[self.integration_volume]
        ]
        return int(vol * (n + 1) * (n + 2) / 2)
