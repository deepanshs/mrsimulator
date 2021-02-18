# -*- coding: utf-8 -*-
"""Base ConfigSimulator class."""
# from mrsimulator.sandbox import AveragingScheme
from pydantic import BaseModel
from pydantic import Field
from typing_extensions import Literal

__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"

# decompose spectrum
__decompose_spectrum_enum__ = {"none": 0, "spin_system": 1}

# integration volume
__integration_volume_enum__ = {"octant": 0, "hemisphere": 1}
__integration_volume_octants__ = [1, 4]

__integration_type_enum__ = {"Alderman": 0, "spherical": 1}


class ConfigSimulator(BaseModel):
    r"""
    The configurable attributes for the Simulator class used in simulation.

    Attributes
    ----------

    number_of_sidebands: int (optional).
        The value is the requested number of sidebands that will be computed in the
        simulation. The value cannot be zero or negative. The default value
        is 64.

    integration_type: enum (optional).
        Orientation sampling type. The valid literals of this enumeration are

        - ``Spherical`` (default), and
        - ``Aldermans``

    integration_volume: enum (optional).
        The volume over which the solid-state spectral frequency integral is evaluated.
        The valid literals of this enumeration are

        - ``octant`` (default), and
        - ``hemisphere``

    integration_density: int (optional).
        The integration density or equivalently the number of orientations sampled over
        an octant. If :math:`n` is the integration_density, then the total number of
        orientation, :math:`N`, sampled over an octant is

        .. math::
            N = \frac{(n+1)(n+2)}{2}.

        The default value is 70.

    decompose_spectrum: enum (optional).
        The value specifies how a simulation result is decomposed into an array of
        spectra. The valid literals of this enumeration are

        - ``none`` (default): When the value is `none`, the resulting simulation is a
          single spectrum, which is an integration of the spectra over all spin systems.
        - ``spin_system``:  When the value is `spin_system`, the resulting simulation
          is an array of spectra, where each spectrum arises from a spin system within
          the Simulator object.

    Example
    -------

    >>> a = Simulator()
    >>> a.config.number_of_sidebands = 128
    >>> a.config.integration_type = 'Alderman'
    >>> a.config.integration_density = 96
    >>> a.config.integration_volume = 'hemisphere'
    >>> a.config.decompose_spectrum = 'spin_system'
    """

    number_of_sidebands: int = Field(default=64, gt=0)
    integration_type: Literal["Alderman", "spherical"] = "Alderman"
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
        py_dict["integration_type"] = __integration_type_enum__[self.integration_type]
        return py_dict

    # averaging scheme. This contains the c pointer used in frequency evaluation
    # at c-level.
    # @property
    # def averaging_scheme(self):
    #     return self._averaging_scheme

    # @averaging_scheme.setter
    # def averaging_scheme(self, other):
    #     if isinstance(other, AveragingScheme):
    #         self._averaging_scheme = other
    #     raise ValueError("Expecting an instance of either the AveragingScheme class.")

    def get_orientations_count(self):
        """Return the total number of orientations.

        Example
        -------

        >>> a = Simulator()
        >>> a.config.integration_density = 20
        >>> a.config.integration_volume = 'hemisphere'
        >>> a.config.get_orientations_count() # (4 * 21 * 22 / 2) = 924
        924
        """
        n = self.integration_density
        vol = __integration_volume_octants__[
            __integration_volume_enum__[self.integration_volume]
        ]
        return int(vol * (n + 1) * (n + 2) / 2)
