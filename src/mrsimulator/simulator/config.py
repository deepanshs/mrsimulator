"""Base ConfigSimulator class."""
# from mrsimulator.sandbox import AveragingScheme
from enum import Enum
from typing import Optional

import numpy as np
from mrsimulator.utils.parseable import Parseable
from pydantic.v1 import BaseModel
from pydantic.v1 import Field
from typing_extensions import Literal


__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"

# decompose spectrum
__decompose_spectrum_enum__ = {"none": 0, "spin_system": 1}
__isotropic_interpolation_enum__ = {"linear": 0, "gaussian": 1}

# integration volume
__integration_volume_enum__ = {"octant": 0, "hemisphere": 1, "sphere": 2}
__integration_volume_octants__ = [1, 4, 8]


class StrType(str, Enum):
    simpson: str = "simpson"
    default: str = "default"


class CustomSampling(BaseModel):
    r"""Custom orientational sampling for powder averaging.

    Attributes
    ----------

    alpha: ndarray
        An array of size N of :math:`\alpha` angle coordinates in radians.

    beta: ndarray
        An array of size N of :math:`\beta` angle coordinates in radians.

    weight: ndarray
        An array of size N of weights corresponding to :math:`(\alpha, \beta)`
        coordinates.

    vertex_indexes: ndarray (optional)
        A 2D array of shape (N, 3) with each row listing the three indexes of
        :math:`(\alpha, \beta)` coordinates that forms a triangular mesh on a unit
        sphere. The indexes are used in frequency interpolation. The default value
        is None and corresponds to no interpolation.
    """
    alpha: Optional[np.ndarray] = None
    beta: Optional[np.ndarray] = None
    weight: Optional[np.ndarray] = None
    vertex_indexes: Optional[np.ndarray] = None

    class Config:
        extra = "forbid"
        allow_population_by_field_name = True
        validate_assignment = True
        arbitrary_types_allowed = True

    def save(
        self, filename: str, target: StrType = StrType.default, units: str = "rad"
    ):
        if units == "rad":
            array = np.array([self.alpha, self.beta, self.weight])
        if units == "deg":
            fn = np.rad2deg
            array = np.array([fn(self.alpha), fn(self.beta), self.weight])
        header = (
            str(array.shape[1]) if target == StrType.simpson else "alpha beta gamma"
        )
        np.savetxt(filename, array.T, header=header)


class ConfigSimulator(Parseable):
    r"""The configurable attributes for the Simulator class used in the simulation.

    Attributes
    ----------

    number_of_sidebands: int (optional).
        Number of sidebands to evaluate in the simulation. The default value is 64.
        Value cannot be negative or zero.

    number_of_gamma_angles: int (optional).
        Number of gamma angles averages in the simulation. The default value is 1. Value
        cannot be negative or zero.

    integration_volume: enum (optional).
        The spatial volume over which the spectral frequency integration/averaging
        is performed. The valid literals of this enumeration are

        - ``octant`` (default),
        - ``hemisphere``, and
        - ``sphere``

    integration_density: int (optional).
        The integration/sampling density or equivalently the number of (alpha, beta)
        orientations over which the frequency spatial averaging is performed within the
        given volume. If :math:`n` is the integration_density, then the total number of
        orientations is given as

        .. math::
            n_\text{octants} \frac{(n+1)(n+2)}{2} n_\gamma,

        where :math:`n_\text{octants}` is the number of octants in the given volume and
        :math:`n_\gamma` is the number of gamma angles. The default value is 70.

    decompose_spectrum: enum (optional).
        The value specifies how a simulation result is decomposed into an array of
        spectra. The valid literals of this enumeration are

        - ``none`` (default): When the value is `none`, the resulting simulation is a
          single spectrum, which is an integration of the spectra over all spin systems.
        - ``spin_system``:  When the value is `spin_system`, the resulting simulation
          is an array of spectra, where each spectrum arises from a spin system within
          the Simulator object.

    isotropic_interpolation: enum (optional).
        Interpolation scheme for isotropic binning. The valid literals are

        - ``linear`` (default): linear interpolation.
        - ``gaussian``:  Gaussian interpolation with `sigma=0.25*bin_width`.

    custom_sampling: CustomSampling
        A CustomSampling object specifying the coordinates and weights used in powder
        averaging.

    is_complex: bool
        When True, the simulation computes both real and imaginary components of the
        signal. When False, only the real component is computed. Default is True.

    Example
    -------

    >>> a = Simulator()
    >>> a.config.number_of_sidebands = 128
    >>> a.config.number_of_gamma_angles = 10
    >>> a.config.integration_density = 96
    >>> a.config.integration_volume = 'hemisphere'
    >>> a.config.decompose_spectrum = 'spin_system'
    """

    number_of_sidebands: int = Field(default=64, gt=0)
    number_of_gamma_angles: int = Field(default=1, gt=0)
    integration_volume: Literal["octant", "hemisphere", "sphere"] = "octant"
    integration_density: int = Field(default=70, gt=0)
    decompose_spectrum: Literal["none", "spin_system"] = "none"
    isotropic_interpolation: Literal["linear", "gaussian"] = "linear"
    custom_sampling: Optional[CustomSampling] = None
    is_complex: bool = True

    class Config:
        extra = "forbid"
        allow_population_by_field_name = True
        validate_assignment = True
        arbitrary_types_allowed = True

    def get_int_dict(self):
        py_dict = self.dict(
            exclude={
                "property_units",
                "name",
                "description",
                "label",
                "custom_sampling",
            }
        )
        py_dict["integration_volume"] = __integration_volume_enum__[
            self.integration_volume
        ]
        py_dict["decompose_spectrum"] = __decompose_spectrum_enum__[
            self.decompose_spectrum
        ]
        py_dict["isotropic_interpolation"] = __isotropic_interpolation_enum__[
            self.isotropic_interpolation
        ]
        if self.custom_sampling is not None:
            py_dict["alpha"] = self.custom_sampling.alpha
            py_dict["beta"] = self.custom_sampling.beta
            py_dict["weight"] = self.custom_sampling.weight
            if self.custom_sampling.vertex_indexes is not None:
                py_dict["positions"] = self.custom_sampling.vertex_indexes.ravel()
            else:
                py_dict["interpolation"] = False
            py_dict["user_defined"] = True
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
        if self.custom_sampling is not None:
            return self.number_of_gamma_angles * self.custom_sampling.alpha.size
        n = self.integration_density
        vol = __integration_volume_octants__[
            __integration_volume_enum__[self.integration_volume]
        ]
        return int(vol * (n + 1) * (n + 2) / 2) * self.number_of_gamma_angles
