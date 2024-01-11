"""Base ConfigSimulator class."""
# from mrsimulator.sandbox import AveragingScheme
from mrsimulator.utils.parseable import Parseable
from pydantic import Field
from typing_extensions import Literal

__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"

# decompose spectrum
__decompose_spectrum_enum__ = {"none": 0, "spin_system": 1}
__isotropic_interpolation_enum__ = {"linear": 0, "gaussian": 1}

# integration volume
__integration_volume_enum__ = {"octant": 0, "hemisphere": 1, "sphere": 2}
__integration_volume_octants__ = [1, 4, 8]


class ConfigSimulator(Parseable):
    r"""The configurable attributes for the Simulator class used in simulation.

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
        orientation is given as

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

    class Config:
        extra = "forbid"
        allow_population_by_field_name = True
        validate_assignment = True

    def get_int_dict(self):
        py_dict = self.dict(exclude={"property_units", "name", "description", "label"})
        py_dict["integration_volume"] = __integration_volume_enum__[
            self.integration_volume
        ]
        py_dict["decompose_spectrum"] = __decompose_spectrum_enum__[
            self.decompose_spectrum
        ]
        py_dict["isotropic_interpolation"] = __isotropic_interpolation_enum__[
            self.isotropic_interpolation
        ]
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
        return int(vol * (n + 1) * (n + 2) / 2) * self.number_of_gamma_angles
