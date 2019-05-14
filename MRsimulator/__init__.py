

from .utils import __get_spin_attribute__
import numpy as np
from .unit import string_to_quantity

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.90@osu.edu"


class Spectrometer:
    """
    Parameters
    ----------
    magnetic_flux_density: str
        The magnetic flux density of the spectromenter.
    sample_rotation_frequency: str
        The rotation frequency of the sample holder.
    sample_rotation_axis: dict
        A dictionary describing the rotation axis of the sample holder with
        keys, `polar_angle` and `azimuthal_angle`. The polar angle is angle
        of the sample holder with respect to the lab frame z-axis while the
        `azimuthal_angle` is the angle with respect to the lab frame x-axis.

    .. doctest::

        >>> from mrsimulator import spectrometer
        >>> s1 = spectrometer(
        ...             magnetic_flux_density='9.4 T',
        ...             sample_rotation_frequency='10 kHz',
        ...             sample_rotation_axis={
        ...                 'polar_angle': '54.735 deg',
        ...                 'azimuthal_angle': '0 deg',
        ...             }
        ...       )
    
    """

    __slots__ = (
        'magnetic_flux_density',
        'sample_rotation_frequency',
        'sample_rotation_axis',
        '_sample',
        '_detector'
    )

    def __init__(self,
                 magnetic_flux_density,
                 sample_rotation_frequency,
                 sample_rotation_axis):
        """Initialize"""

        self.magnetic_flux_density = string_to_quantity(
                                        magnetic_flux_density
                                    ).to('T')
        self.sample_rotation_frequency = string_to_quantity(
                                        sample_rotation_frequency
                                    )
        self.sample_rotation_axis = {
            'polar_angle': string_to_quantity(
                sample_rotation_axis['polar_angle']
            ).to('rad'),
            'azimuthal_angle': string_to_quantity(
                sample_rotation_axis['azimuthal_angle']
            ).to('rad')
        }

    def __value(self):
        return {
            'magnetic_flux_density': self.magnetic_flux_density,
            'sample_rotation_frequency': self.sample_rotation_frequency,
            'sample_rotation_axis': self.sample_rotation_axis
        }

    def load_sample(self, sample):
        self._sample = System(**sample)

    def set_detector(self, detector):
        self._detector = detector

    def run(self, method, sample=None, detector=None, *args, **kwargs):
        if sample is not None:
            sample = System(**sample)
        else:
            sample = self._sample
        
        if detector is None:
            detector = self._detector

        return method(spectrometer=self.__value(), 
                      sites=sample,
                      observed=detector, *args, **kwargs)


class System:
    __slots__ = ()

    def __new__(self,
                sites=[],
                couplings=[]):
        if not isinstance(sites, list):
            raise ValueError((
                f"Expecting a list of sites. Found {type(sites)}."
            ))
        _sites = []
        for site in sites:
            _sites.append(Site(**site))
        return _sites

class Site:
    __slots__ = ()

    def __new__(self,
                isotope_symbol='1H',
                isotropic_chemical_shift='0 Hz',
                shielding_symmetric=None):
        """Initialize."""
        return {
            'isotope_symbol': isotope_symbol,
            'isotropic_chemical_shift': string_to_quantity(
                                            isotropic_chemical_shift
                                        ).to('Hz').value,
            'shielding_symmetric': {
                'anisotropy': string_to_quantity(
                                shielding_symmetric['anisotropy']
                            ).to('Hz').value,
                'asymmetry': float(shielding_symmetric['asymmetry'])
            }
        }
        