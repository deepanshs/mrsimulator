
from .parameterized_tensor import (
    NuclearShieldingTensor,
    ElectricQuadrupoleTensor
)
import json
from .core import (
    _get_json_object,
    _initialize_class_attributes
)
from .utils import __get_spin_attribute__
from .unit import string_to_quantity
from astropy import units as u
from astropy.units import Quantity
from copy import deepcopy


__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.90@osu.edu"


class SpinAveraging:
    __slots__ = (
        '_spin_frequency',
        '_rotor_angle',
        '_magnetic_field',
        '_temperature',
        '_orientation_type',
        '_orientation_'
    )

    def __init__(self,
                 spin_frequency = '0 Hz', 
                 rotor_angle = '54.7356 deg',
                 magnetic_field = '9.4 T',
                 orientation = 'triangle_interpolation(32)'):
        r"""Initialize and return an instance of HaeberlenConvension class."""

        _dictionary_spin_average = {
            'spin_frequency': 0*u.Hz,
            'rotor_angle': 54.7356*u.deg,
            'magnetic_field': 9.4*u.T
        }

        _dictionary_spin_average = _initialize_class_attributes(
            _dictionary_spin_average, spin_frequency,
            rotor_angle, magnetic_field
        )

        m_ = _dictionary_spin_average['magnetic_field']
        if m_.unit.physical_type != 'magnetic flux density':
            raise Exception(f"{m_} is an invalid value for magnetic field.")
        self._magnetic_field = m_

        m_ = _dictionary_spin_average['rotor_angle']
        if m_.unit.physical_type != 'angle':
            raise Exception(f"{m_} is an invalid value for rotor angle.")
        self._rotor_angle = m_

        m_ = _dictionary_spin_average['spin_frequency']
        if m_.unit.physical_type != 'frequency':
            raise Exception(f"{m_} is an invalid value for spin_frequency.")
        self._spin_frequency = m_

        # if 



class System:
    __slots__ = (
        '_sites',
        '_couplings',
        '_spin_frequency',
        '_rotor_angle',
        '_magnetic_field',
        '_temperature',
        '_orientation_type',
        '_orientation_values'
    )

    __public__ = (
        'spin_frequency',
        'rotor_angle',
        'magnetic_field',
        'sites',
    )
    def __init__(self,
                 magnetic_field = '9.4 T',
                 spin_frequency = '0 Hz',
                 rotor_angle = '54.7356 deg',
                 sites = [],
                 couplings = []
                 ):
        r"""Initialize and return an instance of HaeberlenConvension class."""

        _dictionary_system = {
            'spin_frequency': 0*u.Hz,
            'rotor_angle': 54.7356*u.deg,
            'magnetic_field': 9.4*u.T
        }

        self._sites = []
        _dictionary_system = _initialize_class_attributes(
            _dictionary_system, spin_frequency,
            rotor_angle, magnetic_field
        )

        m_ = _dictionary_system['magnetic_field']
        if m_.unit.physical_type != 'magnetic flux density':
            raise Exception(f"{m_} is an invalid value for magnetic field.")
        self._magnetic_field = m_

        m_ = _dictionary_system['rotor_angle']
        if m_.unit.physical_type != 'angle':
            raise Exception(f"{m_} is an invalid value for rotor angle.")
        self._rotor_angle = m_

        m_ = _dictionary_system['spin_frequency']
        if m_.unit.physical_type != 'frequency':
            raise Exception(f"{m_} is an invalid value for spin_frequency.")
        self._spin_frequency = m_

        for i in range(len(sites)):
            self._sites.append( Spin(**sites[i]) )
            self._sites[-1]._larmor_frequency = self._sites[-1]._gyromagnetic_ratio*self._magnetic_field



    @property
    def sites(self):
        """Return a list of Spin instances."""
        return self._sites

# magnetic field
    @property
    def magnetic_field(self):
        """The external magnetic field, :math:`B_z`."""
        return self._magnetic_field
    
    @magnetic_field.setter
    def magnetic_field(self, value):
        if isinstance(value, str):
            value = string_to_quantity(value)
        if self._magnetic_field.unit.physical_type == value.unit.physical_type:
            self._magnetic_field = value
            for i in range(len(self._sites)):
                self._sites[i]._larmor_frequency = self._sites[i]._gyromagnetic_ratio*self._magnetic_field

# spinning frequency
    @property
    def spin_frequency(self):
        """The sampling spinning frequency."""
        return self._spin_frequency
    
    @spin_frequency.setter
    def spin_frequency(self, value):
        if isinstance(value, str):
            value = string_to_quantity(value)
        if self._spin_frequency.unit.physical_type == value.unit.physical_type:
            self._spin_frequency = value

# rotor angle
    @property
    def rotor_angle(self):
        """The sample spin axis with respect to the z-axis."""
        return self._rotor_angle
    
    @rotor_angle.setter
    def rotor_angle(self, value):
        if isinstance(value, str):
            value = string_to_quantity(value)
        if self._rotor_angle.unit.physical_type == value.unit.physical_type:
            self._rotor_angle = value

# methods
    def __repr__(self):
        return f"<System {str(self._get_dict_())}>"
    
    def add_sites(self, value):
        self._sites.append(Spin(value))

    def _get_dict_(self):
        d_ = {}
        for key in self.__class__.__public__:
            value = self.__getattribute__(key)
            if key in ['sites']:
                d_[key] = []
                for i in range(len(value)):
                    if type(value[i]) is Spin:
                        d_[key].append(value[i]._get_dict_())
            elif isinstance(value, Quantity):
                d_[key] = str(value)
            elif value is not ():
                d_[key] = value
        return d_

    def dumps(self):
        """Return the instance as a json serialized string."""

        # d_ = {}
        # for key in self.__class__.__public__:
        #     value = self.__getattribute__(key)
        #     if key in ['sites']:
        #         d_[key] = []
        #         for i in range(len(value)):
        #             if type(value[i]) is Spin:
        #                 d_[key].append(value[i]._get_dict_())
        #     elif value is not ():
        #         d_[key] = value

        return json.dumps(self._get_dict_(), indent=4, allow_nan=False, ensure_ascii=True)


class Spin:
    """
    The class creates an instance of a nuclear spin with the given
    spin isotope symbol.

    For example, an instance of '1H' nucleus may be initialized using

    .. doctest::
    
        >>> from mrsimulator.spin import Spin
        >>> H1 = Spin(isotope_symbol = '1H',
        ...           nuclear_shielding_tensor={
        ...               'haeberlen':('1 Hz', '10 kHz', 0)
        ...           }
        ...      )
        >>> print(H1.dumps())
        {
            "isotope_symbol": "1H",
            "nuclear_shielding_tensor": {
                "haeberlen": {
                    "iso": "1.0 Hz",
                    "zeta": "10.0 kHz",
                    "eta": 0.0
                },
                "euler_angles": {
                    "alpha": "0.0 rad",
                    "beta": "0.0 rad",
                    "gamma": "0.0 rad"
                }
            }
        }

    Likewise, an instance of '23Na' nucleus may be initialized using

    .. doctest::

        >>> Na23 = Spin(isotope_symbol = '23Na', 
        ...             nuclear_shielding_tensor={
        ...                 'principal_components': ('12 kHz', '10 kHz', '-20 kHz'),
        ...                 'euler_angles': ('0 deg', '54.74 deg')
        ...             },
        ...             electric_quadrupole_tensor={
        ...                 'Cq': '1.2 MHz',
        ...                 'eta': 0.1
        ...             }
        ...         )
        >>> print(Na23.dumps())
        {
            "isotope_symbol": "23Na",
            "nuclear_shielding_tensor": {
                "haeberlen": {
                    "iso": "0.6666666666666666 kHz",
                    "zeta": "-20.666666666666668 kHz",
                    "eta": 0.0967741935483871
                },
                "euler_angles": {
                    "alpha": "0.0 deg",
                    "beta": "54.74 deg",
                    "gamma": "0.0 rad"
                }
            },
            "electric_quadrupole_tensor": {
                "Cq": "1.2 MHz",
                "eta": 0.1,
                "euler_angles": {
                    "alpha": "0.0 rad",
                    "beta": "0.0 rad",
                    "gamma": "0.0 rad"
                }
            }
        }

    .. note::
        The value of arguments ``nuclear_shielding_tensor`` and
        ``electric_quadrupole_tensor`` are forwarded to
        :ref:`shielding_tensor_api` and :ref:`quad_tensor_api` class
        respectively.
    """

    __slots__= (
        '_isotope_symbol',
        '_spin_quantum_number',
        '_gyromagnetic_ratio',
        '_nuclear_shielding_tensor',
        '_larmor_frequency',
        '_electric_quadrupole_tensor'
    )

    __public__ = (
        'isotope_symbol',
        'nuclear_shielding_tensor',
        'electric_quadrupole_tensor'
    )

    def __init__(
            self, 
            isotope_symbol='1H',
            nuclear_shielding_tensor = (),
            electric_quadrupole_tensor = (),
            *args,
            **kwargs):

        
        self._isotope_symbol = isotope_symbol
        sp_ = __get_spin_attribute__[isotope_symbol]
        self._spin_quantum_number = sp_['spin']
        self._gyromagnetic_ratio = string_to_quantity(str(sp_['gyromagnetic_ratio']) + ' MHz/T')
        self._larmor_frequency = 50000*u.nT*self._gyromagnetic_ratio

        if type(nuclear_shielding_tensor) is NuclearShieldingTensor:
            self._nuclear_shielding_tensor = nuclear_shielding_tensor
        else:
            self._nuclear_shielding_tensor = NuclearShieldingTensor(nuclear_shielding_tensor)
        
        self._electric_quadrupole_tensor = ()
        if self._spin_quantum_number > 0.5:
            if type(electric_quadrupole_tensor) is ElectricQuadrupoleTensor:
                self._electric_quadrupole_tensor = electric_quadrupole_tensor
            else:
                self._electric_quadrupole_tensor = ElectricQuadrupoleTensor(electric_quadrupole_tensor)
            
# larmor frequency
    @property
    def larmor_frequency(self):
        r"""The nuclear larmor frequency, :math:`\nu = -\frac{1}{2\pi}\gamme_n B_z`."""
        return self._larmor_frequency

        
# isotope symbol
    @property
    def isotope_symbol(self):
        """The isotope.
        
        .. doctest::

            >>> H1.isotope_symbol
            '1H'

            >>> Na23.isotope_symbol
            '23Na'

        """
        return self._isotope_symbol

# spin quantum number
    @property
    def spin_quantum_number(self):
        r"""The spin quantum number, :math:`s`, asoociated with the isotope.
        
        .. doctest::

            >>> H1.spin_quantum_number
            0.5

            >>> Na23.spin_quantum_number
            1.5

        """
        return self._spin_quantum_number
    
# gyromagnetic ratio
    @property
    def gyromagnetic_ratio(self):
        r"""The gyromagnetic ratio of the nucleus, :math:`\gamma_n`, associated with the nuclear spin.
        
        .. doctest::

            >>> H1.gyromagnetic_ratio
            <Quantity 42.57747892 MHz / T>

            >>> Na23.gyromagnetic_ratio
            <Quantity 11.262 MHz / T>
        
        """
        return self._gyromagnetic_ratio
    
    @property
    def nuclear_shielding_tensor(self):
        r"""
        Return the nuclear shielding tensor of the nucleus as a
        :ref:`shielding_tensor_api` instance.
        """
        return self._nuclear_shielding_tensor
    
    @property
    def electric_quadrupole_tensor(self):
        r"""
        Return the electric quadrupole tensor of the nucleus as a
        :ref:`quad_tensor_api` instance.
        """
        return self._electric_quadrupole_tensor

    def _get_dict_(self):
        d_ = {}
        for key in self.__class__.__public__:
            value = self.__getattribute__(key)
            if type(value) in [NuclearShieldingTensor, ElectricQuadrupoleTensor]:
                d_[key] = _get_json_object(value.__dict_obj__())
            elif value is not ():
                d_[key] = value
        return d_

    def dumps(self):
        """
        Return the instance as a json serialized string.
        
        .. doctest::

            >>> print(H1.dumps())
            {
                "isotope_symbol": "1H",
                "nuclear_shielding_tensor": {
                    "haeberlen": {
                        "iso": "1.0 Hz",
                        "zeta": "10.0 kHz",
                        "eta": 0.0
                    },
                    "euler_angles": {
                        "alpha": "0.0 rad",
                        "beta": "0.0 rad",
                        "gamma": "0.0 rad"
                    }
                }
            }

        """
        # d_ = {}
        # for key in self.__class__.__public__:
        #     value = self.__getattribute__(key)
        #     if type(value) in [NuclearShieldingTensor, ElectricQuadrupoleTensor]:
        #         d_[key] = _get_json_object(value.__dict_obj__())
        #     elif value is not ():
        #         d_[key] = value
                                           
        return json.dumps(self._get_dict_(), indent=4, allow_nan=False, ensure_ascii=True)