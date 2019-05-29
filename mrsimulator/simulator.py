
from .utils import __get_spin_attribute__
from .unit import string_to_quantity, _ppm
import json
from urllib.parse import urlparse
from ._utils_download_file import _download_file_from_url
from . import examples
from astropy import units as u
from copy import deepcopy
try:
    import csdfpy as cp
except ImportError:
    pass


__author__ = "Deepansh J. Srivastava"
__email__ = ["srivastava.89@osu.edu", "deepansh2012@gmail.com"]


def _fn_(x):
    return int(''.join([i for i in x if i.isnumeric()]))


def _import_json(filename):
    res = urlparse(filename)
    if res[0] not in ['file', '']:
        filename = _download_file_from_url(filename)
    with open(filename, "rb") as f:
        content = f.read()
        return (json.loads(str(content, encoding="UTF-8")))


class _Dimensions:
    __slots__ = ()

    def __new__(
            number_of_points=1024,
            spectral_width='100 kHz',
            reference_offset='0 Hz'):

        """Initialize."""
        dictionary = {}
        dictionary['number_of_points'] = int(number_of_points)

        spectral_width = string_to_quantity(spectral_width)
        if spectral_width.unit.physical_type != 'frequency':
            raise Exception((
                "A frequency value is required for the 'spectral_width'."
            ))
        dictionary['spectral_width'] = spectral_width.to('Hz').value

        reference_offset = string_to_quantity(reference_offset)
        if reference_offset.unit.physical_type != 'frequency':
            raise Exception((
                "A frequency value is required for the 'reference_offset'."
            ))
        dictionary['reference_offset'] = reference_offset.to('Hz').value
        return dictionary


class _Spectrum(_Dimensions):
    """Set up a virtual spin environment."""
    __slots__ = ()

    def __new__(
            self,
            magnetic_flux_density='9.4 T',
            rotor_frequency='0 kHz',
            rotor_angle='54.735 deg',
            rotor_phase='0 rad',
            nucleus='1H',
            *args, **kwargs):
        """Initialize"""
        dimension_dictionary = super(_Spectrum, self).__new__(*args, **kwargs)
        magnetic_flux_density = string_to_quantity(magnetic_flux_density)
        if magnetic_flux_density.unit.physical_type != 'magnetic flux density':
            raise Exception((
                "A magnetic flux density quantity is required for "
                "'magnetic_flux_density'."
            ))
        magnetic_flux_density = magnetic_flux_density.to('T').value

        rotor_frequency = string_to_quantity(rotor_frequency)
        if rotor_frequency.unit.physical_type != 'frequency':
            raise Exception((
                "A frequency quantity is required for 'rotor_frequency'."
            ))
        rotor_frequency = rotor_frequency.to('Hz').value

        rotor_angle = string_to_quantity(rotor_angle).to('rad').value
        rotor_phase = string_to_quantity(rotor_phase).to('rad').value

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
        return dictionary


class Isotopomers(list):
    def __init__(self, isotopomers: list) -> list:
        isotopomers_ = _Isotopomers(isotopomers)
        list.__init__(self, isotopomers_)

    def append(self, value):
        isotopomer_ = _Isotopomer(**value)
        list.append(isotopomer_)


class _Isotopomers:
    __slots__ = ()

    def __new__(self, isotopomers: list) -> list:

        if not isinstance(isotopomers, list):
            raise Exception((
                f"A list of isotopomers is required, "
                f"found {type(isotopomers)}."
            ))
        if len(isotopomers) != 0:
            if not isinstance(isotopomers[0], dict):
                raise Exception((
                    f"A list of isotopomer dictionaries is "
                    f"required, found {type(isotopomers[0])}."
                ))

        isotopomers_ = []
        for isotopomer in isotopomers:
            isotopomer_set = _Isotopomer(**isotopomer)
            isotopomers_.append(isotopomer_set)

        return isotopomers_


class _Isotopomer:
    __slots__ = ()

    def __new__(
            self,
            sites: list = [],
            couplings: list = [],
            abundance: str = '100.0 %') -> list:
        """Initialize."""
        if not isinstance(sites, list):
            raise ValueError((
                f"Expecting a list of sites. Found {type(sites)}."
            ))
        _sites = []
        abundance = string_to_quantity(abundance).to('').value
        for site in sites:
            _sites.append(_Site(**site))
        return ({'sites': _sites, 'abundance': abundance})


def _check_values_in_ppm(value, property):
    value_ = string_to_quantity(value)
    if value_.unit.physical_type == 'dimensionless':
        if str(value_.unit) == 'ppm':
            return value_
        else:
            return value_.to(_ppm)
    if value_.unit.physical_type == 'frequency':
        return value_
    else:
        raise Exception(
            (f"Expecting '{property}' in units of frequency or a "
             f"dimensionless frequency ratio, ppm, found {str(value_)}.")
        )


class _Site:
    __slots__ = ()

    def __new__(self,
                isotope_symbol='1H',
                isotropic_chemical_shift='0 ppm',
                shielding_symmetric=None):
        """Initialize."""
        return {
            'isotope_symbol': isotope_symbol,
            'isotropic_chemical_shift': _check_values_in_ppm(
                                            isotropic_chemical_shift,
                                            'isotropic_chemical_shift'),
            'shielding_symmetric': {
                'anisotropy': _check_values_in_ppm(
                                shielding_symmetric['anisotropy'],
                                'shielding anisotropy'),
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
    """

    __slots__ = (
        '_isotopomers_c',
        '_isotopomers',
        '_spectrum',
        '_spectrum_c',
        '_isotope_list',
        '_allowed_isotopes',
        '_larmor_frequency',
        '_freq',
        '_amp'
    )

    def __init__(self, isotopomers=None, spectrum=None):
        self._isotopomers_c = []
        self._isotopomers = []
        self._spectrum = {}
        self._spectrum_c = {}
        self._isotope_list = []
        self._freq = []*u.Hz
        self._amp = []
        isotope_list = __get_spin_attribute__.keys()
        self._allowed_isotopes = list(set([
            isotope for isotope in isotope_list
            if __get_spin_attribute__[isotope]['spin'] == 0.5
        ]))

        if isotopomers is not None:
            self.isotopomers = isotopomers

        if spectrum is not None:
            self.spectrum = spectrum

    @property
    def isotope_list(self):
        """
        Return a list of unique isotopes symbols from the list of
        isotopomers.
        """
        return self._isotope_list

    @property
    def isotopomers(self):
        """
        Return a list of :ref:`isotopomer` objects. The attribute can also be
        used to assign a list of valid isotopomers.
        """
        # return json.dumps(self._isotopomers, ensure_ascii=True, indent=2)
        return self._isotopomers

    @isotopomers.setter
    def isotopomers(self, value):
        self._isotopomers_c = _Isotopomers(value)
        isotope_list = [
            site['isotope_symbol'] for isotopomer in self._isotopomers_c
            for site in isotopomer['sites']
            if site['isotope_symbol'] in self._allowed_isotopes
        ]

        self._isotope_list = list(set(isotope_list))
        self._isotope_list.sort(key=_fn_)
        self._isotopomers = value

    @property
    def spectrum(self):
        """
        Return a :ref:`spectrum` object. The attribute can also be
        used to assign a valid spectrum object.
        """
        # return json.dumps(self._spectrum, ensure_ascii=True, indent=2)
        return self._spectrum

    @spectrum.setter
    def spectrum(self, value):
        self._spectrum_c = _Spectrum(**value['direct_dimension'])
        self._spectrum = value

    def run(self, method, data_object=False, **kwargs):
        """
        Simulate the spectrum using the specified method. The keyword argument
        are the arguments of the specified `method`.

        :ivar method: The method used in computing the linshape.
        :ivar data_object: If true, returns a `csdm` data object. If false,
            returns a tuple of frequency array and the
            corresponding amplitude array. The amplitude is a
            `numpy <https://docs.scipy.org/doc/numpy/reference/generated/numpy.array.html>`_
            array, whereas, the frequency is a
            `Quantity <http://docs.astropy.org/en/stable/units/quantity.html>`_
            array. The default value is False.

        :returns: A `csdm` object if `data_object` is True, else a tuple of
            frequency and amplitude. For details, refer to the
            description of `data_object`.
        """
        isotopomers = self._isotopomers_c
        if isotopomers is []:
            raise Exception("Isotopomers are required for simulation.")
        spectrum = self._spectrum_c
        if spectrum is {}:
            raise Exception((
                "Cannot simulate without the spectrum information."
            ))
        self._freq, self._amp, \
            self._larmor_frequency, list_index_isotopomer = method(
                    spectrum=spectrum,
                    isotopomers=isotopomers, **kwargs
                )

        """The frequency is in the units of Hz."""
        self._freq *= u.Hz
        """The larmor_frequency is in the units of MHz."""
        self._larmor_frequency *= u.MHz

        isotopo_ = [isotopomers[i] for i in list_index_isotopomer]
        application = {
            'isotopomers': str(isotopo_),
            'spectrum': spectrum
        }
        data_object = False
        if data_object:
            return get_csdfpy_object(
                        self._freq,
                        self._larmor_frequency,
                        self._amp,
                        application
            )
        else:
            return deepcopy(self._freq), deepcopy(self._amp)

    def load_isotopomers(self, filename):
        """
        Load a JSON serialized isotopomers file.

        See an
        `example <https://raw.githubusercontent.com/DeepanshS/mrsimulator-test/master/isotopomers_ppm.json>`_
        of JSON serialized isotopomers file. For details, refer to the
        :ref:`load_isotopomers` section.
        """
        contents = _import_json(filename)
        self.isotopomers = contents['isotopomers']


def get_csdfpy_object(x, x0, y, application):
    ob1 = cp.new()
    d1 = {
        'type': 'linear',
        'number_of_points': x.size,
        'increment': str(x[1] - x[0]),
        'index_zero_value': str(x[0]),
        'origin_offset': str(x0),
        'quantity': 'frequency'
    }
    ob1.add_dimension(d1)
    s1 = {
        'type': 'internal',
        'numeric_type': 'float64',
        'components': [y],
        'component_labels': ['arbitrary units'],
    }
    ob1.add_dependent_variable(s1)

    ob1.dependent_variables[0].application = {
            'mrsimulator': application
        }
    return ob1

# class _dimension_object:
#     __slots__ = (
#         'type',
#         'number_of_points',
#         'increment',
#         'index_zero_value',
#         'origin_offset',
#         'coordinates'
#     )

#     def __init__(self, vector, origin_offset):
#         self.type= 'linear'
#         self.number_of_points = vector.size
#         self.increment = str(vector[1] - vector[0]) + ' Hz'
#         self.index_zero_value = str(vector[0]) + ' Hz'
#         self.origin_offset = str(origin_offset) + ' MHz'
#         self.coordinates = vector

#     def to_ppm(self):
#         return None


def _simulator(
        spectrum,
        method,
        isotopomers,
        **kwargs):
    """Execute the method and acquire a spectrum.

    :param method: A function describing the pulse sequence.
    :param sample: A python dictionary describing the sample;

    :returns: freq: A numpy array of frequency values.
    :returns: amp:  A numpy array of amplitudes corresponding to
                    the frequencies.

    """
    if isotopomers is None:
        raise Exception("No isotopomer found.")
    isotopomers = _Isotopomers(isotopomers)

    spectrum = _Spectrum(**spectrum['direct_dimension'])
    # print(spectrum)

    freq, amp = method(
                    spectrum=spectrum,
                    isotopomers=isotopomers, **kwargs
                )

    return (freq, amp)


def run_test():
    from mrsimulator.methods import one_d_spectrum
    import matplotlib.pyplot as plt
    s1 = Simulator()
    # test 1
    s1.isotopomers, s1.spectrum = examples.csa_static()
    freq, amp = s1.run(one_d_spectrum, verbose=1)

    ax = plt.subplots(1, 2, figsize=(6, 3))[1]
    ax[0].plot(freq, amp)
    # ax[0].plot(ob1.dimensions[0].coordinates,
    #            ob1.dependent_variables[0].components[0])
    # label_ = ob1.dimensions[0].axis_label
    label_ = f'frequency / {freq.unit}'
    ax[0].set_xlabel(label_)

    # test 2
    s1.isotopomers, s1.spectrum = examples.csa_mas()
    freq, amp = s1.run(one_d_spectrum, verbose=1)
    ax[1].plot(freq, amp)
    # ax[1].plot(ob2.dimensions[0].coordinates,
    #            ob2.dependent_variables[0].components[0])
    # label_ = ob1.dimensions[0].axis_label
    label_ = f'frequency / {freq.unit}'
    ax[1].set_xlabel(label_)
    plt.tight_layout()
    plt.show()
