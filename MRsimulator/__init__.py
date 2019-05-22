

from .utils import __get_spin_attribute__
import numpy as np
from .unit import string_to_quantity
import json

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.90@osu.edu"


class _Dimensions:
    __slots__ = ()

    def __new__(
        number_of_points=1024,
        spectral_width='100 kHz',
        reference_offset='0 Hz'
    ):
        """Initialize."""
        dictionary = {}
        dictionary['number_of_points'] = int(number_of_points)

        spectral_width  = string_to_quantity(spectral_width)
        if spectral_width.unit.physical_type != 'frequency':
            raise Exception("A frequency value is required for the 'spectral_width'.") 
        dictionary['spectral_width'] = spectral_width.to('Hz').value

        reference_offset  = string_to_quantity(reference_offset)
        if reference_offset.unit.physical_type != 'frequency':
            raise Exception("A frequency value is required for the 'reference_offset'.") 
        dictionary['reference_offset'] = reference_offset.to('Hz').value
        return dictionary


class _Spectrum(_Dimensions):
    """Set up a virtual spin environment."""
    __slots__ = ()

    def __new__(self,
                 magnetic_flux_density='9.4 T',
                 rotor_frequency='2 kHz',
                 rotor_angle='54.735 deg',
                 rotor_phase='0 rad',
                 nucleus='1H',
                 *args, **kwargs):
        """Initialize"""
        dimension_dictionary = super(_Spectrum, self).__new__(*args, **kwargs)
        # print(dimension_dictionary)
        magnetic_flux_density = string_to_quantity(magnetic_flux_density)
        if magnetic_flux_density.unit.physical_type != 'magnetic flux density':
            raise Exception("A magnetic flux density quantity is required for 'magnetic_flux_density'.")
        magnetic_flux_density = magnetic_flux_density.to('T').value

        rotor_frequency = string_to_quantity(rotor_frequency)
        if rotor_frequency.unit.physical_type != 'frequency':
            raise Exception("A frequency quantity is required for 'rotor_frequency'.")
        rotor_frequency = rotor_frequency.to('Hz').value


        rotor_angle = string_to_quantity(rotor_angle).to('rad').value
        rotor_phase = string_to_quantity(rotor_phase).to('rad').value


        # sample_rotation_axis = {
        #     'polar_angle': string_to_quantity(
        #         sample_rotation_axis['polar_angle']
        #     ).to('rad').value,
        #     'azimuthal_angle': string_to_quantity(
        #         sample_rotation_axis['azimuthal_angle']
        #     ).to('rad').value
        # }

        
        dictionary = {
            'magnetic_flux_density': magnetic_flux_density,
            'rotor_frequency': rotor_frequency,
            'rotor_angle': rotor_angle,
            'rotor_phase': rotor_phase
        }

        dictionary.update(dimension_dictionary)

        detect = get_proper_detector_nucleus(nucleus)
        try:
            spin_dictionary = __get_spin_attribute__[detect]
            spin_dictionary['isotope'] = detect
        except KeyError:
            raise Exception(f"Failed to simulates the {detect} spectrum.")

        dictionary.update(spin_dictionary) 

        # print (dictionary)
        return dictionary



class Isotopomers(list):
    def __init__(self, isotopomers:list) -> list:
        isotopomers_, isotope_list = _Isotopomers(isotopomers)
        list.__init__(self, isotopomers_)

    def append(self, value):
        isotopomer_, isotope_list_ = _Isotopomer(**value)
        list.append(isotopomer_)
        # return (isotopomer_set, isotope_list_)


class _Isotopomers:
    __slots__ = ()

    def __new__(self, isotopomers : list) -> list:
        
        if not isinstance(isotopomers, list):
            raise Exception(f"A list of isotopomers is required, found {type(isotopomers)}.")
        if len(isotopomers) != 0:
            if not isinstance(isotopomers[0], dict):
                raise Exception((
                        f"A list of isotopomer dictionaries is "
                        f"required, found {type(isotopomers[0])}."))

        isotopomers_ = []
        isotope_list = []
        for isotopomer in isotopomers:
            isotopomer_set, isotope_list_ = _Isotopomer(**isotopomer)
            isotopomers_.append(isotopomer_set)

            isotope_list.append(isotope_list_)

        # list.__init__(self, isotopomers_)
        return (isotopomers_, isotope_list)

   


class _Isotopomer:
    __slots__ = ()

    def __new__(self,
                sites : list = [],
                couplings : list = [],
                abundance : float = 1.0) -> list:
        # print('\n', sites)
        if not isinstance(sites, list):
            raise ValueError((
                f"Expecting a list of sites. Found {type(sites)}."
            ))
        _sites = []
        isotope_list = []
        for site in sites:
            # print('\n', site)
            _sites.append(_Site(**site))
            isotope_list.append(_sites[-1]['isotope_symbol'])
        return ({'sites': _sites, 'abundance': abundance}, isotope_list)


class _Site:
    __slots__ = ()

    def __new__(self,
                isotope_symbol='1H',
                # abundance=None,
                isotropic_chemical_shift='0 Hz',
                shielding_symmetric=0):
        """Initialize."""
        # if abundance is None:
        #     abundance = __get_spin_attribute__[isotope_symbol]['natural_abundance']
        # else:
        #     abundance = string_to_quantity(abundance).to('')*100.
        return {
            'isotope_symbol': isotope_symbol,
            # 'abundance': abundance,
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


def get_proper_detector_nucleus(string):
    numeric = '0123456789'
    for i, c in enumerate(string):
        if c in numeric:
            break
    return string[i:]+string[0:i]


class Simulator:
    """
    The simulator class.

    :ivar isotopomers: A list of isotopomer dictionaries.
    :ivar spectrum: A dictionary with metadata concerning the simulation.
    """
    # A Simulator object may be created using
    
    # .. doctest::
        
    #         >>> sim1 = Simulator()
    

    __slots__ = (
        '_isotopomers_list',
        '_isotopomers',
        '_spectrum',
        '_isotope_list'
    )


    def __init__(self, isotopomers=None, spectrum=None):
        self._isotopomers_list = []
        self._isotopomers = []
        self._spectrum = {}
        self._isotope_list = []

        if isotopomers is not None:
            self.isotopomers = isotopomers
        
        if spectrum is not None:
            self.spectrum = spectrum


    @property
    def isotope_list(self):
        """Returns a set of isotopes symbols from the list of isotopomers."""
        return self._isotope_list

    @property
    def isotopomers(self):
        """
        Assign a list of isotopomers to the simulator object.
        For example,

        .. doctest::
            
            >>> isotopomers = [
            ...     {
            ...         "sites": [
            ...             {
            ...                 "isotope_symbol": "13C",
            ...                 "isotropic_chemical_shift": "1 Hz",
            ...                 "shielding_symmetric": {
            ...                     "anisotropy": "-3.89 kHz",
            ...                     "asymmetry": 0.25
            ...                 }
            ...             }
            ...         ]
            ...     }
            ... ]
            >>> sim1.isotopomers = isotopomers
        """
        return self._isotopomers_list

    @isotopomers.setter
    def isotopomers(self, value):
        self._isotopomers, isotope_list = _Isotopomers(value)
        isotope_list = [item for sublist in isotope_list for item in sublist]
        self._isotope_list = set(isotope_list)
        self._isotopomers_list = value
    
    @property
    def spectrum(self):
        """
        Assign a new dictionary of spectrum metadata to the simulator object.
        The argument `value` is a dictionary with the following keywords

        .. doctest::
            
            >>> spectrum = {
            ...     "direct_dimension": {
            ...         "magnetic_flux_density": "9.4 T",
            ...         "rotor_frequency": "5 kHz",
            ...         "rotor_angle": "54.735 deg",
            ...         "rotor_phase": "0 rad"
            ...         "number_of_points": 8192,
            ...         "spectral_width": "100 kHz",
            ...         "reference_offset": "0 Hz",
            ...         "nucleus": "13C"
            ...     }
            ... }
            >>> sim1.spectrum = spectrum
        """
        return self._spectrum

    @spectrum.setter
    def spectrum(self, value):
        self._spectrum = _Spectrum(**value['direct_dimension'])

    def run(self, method, **kwargs):
        """Simulate the spectrum using the specified method."""
        isotopomers = self._isotopomers
        if isotopomers is None:
            raise Exception("Isotopomers are required for simulation.")
        spectrum = self._spectrum
        if spectrum is None:
            raise Exception("Cannot simulate without the spectrum information.")
        freq, amp = method(
                    spectrum=spectrum, 
                    isotopomers=isotopomers, **kwargs
                )
        return (freq, amp)


def _simulator(
        spectrum,
        method,
        isotopomers,
        **kwargs):
    """Execute the method and acquire a spectrum.

    :param method: A function describing the pulse sequence.
    :param sample: A python dictionary describing the sample; 

    :returns: freq: A numpy array of frequency values.
    :returns: amp:  A numpy array of amplitudes corresponding to the frequencies.

    """
    if isotopomers is None:
        raise Exception("No isotopomer found.")
    # if spectrum is None:
    #     raise Exception("No isotopomer found.")
    isotopomers = _Isotopomers(isotopomers)


    spectrum = _Spectrum(**spectrum['direct_dimension'])
    print(spectrum)

    # if 
    # environment = _Environment(**environment)

    # detect = get_proper_detector_nucleus(detect)
    # try:
    #     spin_dictionary = __get_spin_attribute__[detect]
    #     spin_dictionary['isotope'] = detect
    # except KeyError:
    #     raise Exception(f"Failed to simulates the {detect} spectrum.")
    freq, amp = method(
                    spectrum=spectrum, 
                    isotopomers=isotopomers, **kwargs
                )

    return (freq, amp)