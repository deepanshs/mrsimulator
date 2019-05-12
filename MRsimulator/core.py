
import numpy as np
from .unit import string_to_quantity
from astropy.units import Quantity
from astropy import units as u
import json


__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.90@osu.edu"


CLASS_LIST = (
    'str',
    'Quantity',
    'float',
    'float64',
    'int',
    'HaeberlenConvension',
    'PrincipalComponents',
    'EulerAngles',
)


def _get_json_object(dict_, keys=None):
    if keys is None:
        keys = dict_.keys()
    d_ = {}
    for key in keys:
        if isinstance(dict_[key], Quantity):
            d_[key] = str(dict_[key])
        elif type(dict_[key]) in [int, float, str, None]:
            d_[key] = dict_[key]
        else:
            d_[key] = dict_[key].__dict_obj__()
    # d_ = { key: str(dict_[key]) for key in keys }
    return d_


def _initialize_class_attributes(dictionary, *arg, **kwargs):
    if arg == () and kwargs == ():
        return dictionary

    if arg != ():
        if isinstance(arg[0], dict):
            input_dict = arg[0]
        else:
            keys = dictionary.keys()
            input_dict = {key: arg[i] for i, key in enumerate(keys) if i<len(arg)}
    else:
        input_dict = kwargs

    default_keys = dictionary.keys()
    input_keys = input_dict.keys()

    for key in input_keys:
        if key not in default_keys:
            raise KeyError((
                f"Encountered an invalid key, '{key}'. The valid keys are "
                f"{str(list(default_keys))[1:-1]}."
            ))

        if input_dict[key].__class__.__name__ not in CLASS_LIST:
            raise ValueError((
                f"{input_dict[key].__class__.__name__} is not a valid value "
                f"type for {key} key, '{input_dict[key]}'."
            ))

        class_type = type(dictionary[key])

        if not isinstance(input_dict[key], class_type):
            if isinstance(class_type, Quantity):
                dictionary[key] = string_to_quantity(str(input_dict[key]))
            else:
                dictionary[key] = dictionary[key].__class__(input_dict[key])
        else:
            dictionary[key] = input_dict[key]

    return dictionary


def _get_haeberlen_values_from_PAS(dictionary):
    """Convert and return the principal components to HaeberlenConvension."""
    xx = dictionary['xx']
    yy = dictionary['yy']
    zz = dictionary['zz']
    # unit = xx.unit
    iso = (xx + yy + zz)/3.0
    xx = xx - iso
    yy = yy - iso
    zz = zz - iso
    zeta = max(abs(xx), abs(yy), abs(zz)) * np.sign(zz)
    if zeta == 0.:
        eta = 0.
    else:
        if abs(xx) == abs(zeta):
            eta = abs(yy-zz)/abs(zeta)
        elif abs(yy) == abs(zeta):
            eta = abs(xx-zz)/abs(zeta)
        elif abs(zz) == abs(zeta):
            eta = abs(xx-yy)/abs(zeta)
    return iso, zeta, eta


def _get_PAS_from_haeberlen_values(dictionary):
    """Convert and return HaeberlenConvension to principal components."""
    # print('dict', dictionary)
    iso = dictionary['iso']
    zeta = dictionary['zeta']
    eta = dictionary['eta']
    zz = zeta + iso
    temp1 = 0.5*eta*zeta
    temp2 = iso - 0.5*zeta
    yy = temp2 + temp1
    xx = temp2 - temp1
    return xx, yy, zz


class HaeberlenConvension:
    r"""
    An instance of this class represents the components of a second rank
    symmetric tensor in the principal axis frame with three parameters,
    the isotropic :math:`\delta`, anisotropic :math:`\zeta` and the asymmetry
    :math:`\eta` parameters. The
    instance may be initialized using any of the following-

        >>> from mrsimulator.core import HaeberlenConvension, PrincipalComponents, EulerAngles

        >>> a = HaeberlenConvension('0 Hz', '-6 Hz', 0)
        >>> print(a)
        {'iso': <Quantity 0. Hz>, 'zeta': <Quantity -6. Hz>, 'eta': 0.0}

        >>> b = HaeberlenConvension(iso='5 Hz', zeta='-6 Hz', eta=0.6666666666666666)
        >>> print(b)
        {'iso': <Quantity 5. Hz>, 'zeta': <Quantity -6. Hz>, 'eta': 0.6666666666666666}

        >>> c = HaeberlenConvension({'iso': '15 Hz', 'zeta': '1 kHz', 'eta':'0.95'})
        >>> print(c)
        {'iso': <Quantity 15. Hz>, 'zeta': <Quantity 1. kHz>, 'eta': 0.95}

    :returns: An instance of the :ref:`haeberlen_api` class.
    """

    __slots__ = (
        '_dictionary_haeberlen',
    )
    def __init__(self, iso='0 Hz', zeta='0 Hz', eta=0): # *args, **kwargs):
        r"""Initialize and return an instance of HaeberlenConvension class."""
        self._dictionary_haeberlen = {
            'iso': 0*u.Hz,
            'zeta': 0*u.Hz,
            'eta': 0.0
        }

        self._dictionary_haeberlen = _initialize_class_attributes(
            self._dictionary_haeberlen, iso, zeta, eta
        )
        
        l_ = self._dictionary_haeberlen['iso'].unit.physical_type
        r_ = self._dictionary_haeberlen['zeta'].unit.physical_type
        if l_ != r_:
            raise Exception((
                "'iso' and 'zeta' parameters must have the same physical type, found "
                f"{l_} and {r_}, respectively."
            ))

# private methods
    def __repr__(self):
        return f'<HaeberlenConvension {str(self._dictionary_haeberlen)}>'

    def __str__(self):
        return str(self._dictionary_haeberlen)

# attributes
    @property
    def iso(self):
        r"""The isotropic nuclear shielding, :math:`\delta`."""
        return self._dictionary_haeberlen['iso']

    @property
    def zeta(self):
        r"""The nuclear shielding anisotropic, :math:`\zeta`."""
        return self._dictionary_haeberlen['zeta']

    @property
    def eta(self):
        r"""The nuclear shielding asymmetry, :math:`\eta`."""
        return self._dictionary_haeberlen['eta']

# methods
    def principal_values(self):
        r"""
        The principal components corresponding to the Haeberlen values.

        The components are evaluated as

        :math:`zz = \zeta + \delta`

        :math:`yy = \delta - \frac{1}{2}\zeta + \frac{1}{2}\eta\zeta`

        :math:`xx = \delta - \frac{1}{2}\zeta - \frac{1}{2}\eta\zeta`

            >>> print(a.principal_values())
            {'xx': <Quantity 3. Hz>, 'yy': <Quantity 3. Hz>, 'zz': <Quantity -6. Hz>}

            >>> print(b.principal_values())
            {'xx': <Quantity 10. Hz>, 'yy': <Quantity 6. Hz>, 'zz': <Quantity -1. Hz>}

            >>> print(c.principal_values())
            {'xx': <Quantity -960. Hz>, 'yy': <Quantity -10. Hz>, 'zz': <Quantity 1.015 kHz>}

        :return: An instance of :ref:`principal_api` class.
        """
        xx, yy, zz = _get_PAS_from_haeberlen_values(self._dictionary_haeberlen)
        return PrincipalComponents(xx=xx, yy=yy, zz=zz)

    def __dict_obj__(self):
        d_ = _get_json_object(self._dictionary_haeberlen)
        return d_

    def dumps(self):
        """
        Return the instance as a json serialized string.

            >>> print(a.dumps())
            {
                "iso": "0.0 Hz",
                "zeta": "-6.0 Hz",
                "eta": 0.0
            }

        """
        d_ = _get_json_object(self._dictionary_haeberlen)
        return json.dumps(d_, indent=4, allow_nan=False, ensure_ascii=True)

class PrincipalComponents:
    r"""
    An instance of this class represents the principal components of a
    second rank symmetric tensor.

    The principal components are denoted as :math:`xx`, :math:`yy`, and
    :math:`zz` respectrively.
    The instance may be initialized using any of the following-

        >>> a = PrincipalComponents('1 Hz', '15 Hz', '-6 Hz')
        >>> print(a)
        {'xx': <Quantity 1. Hz>, 'yy': <Quantity 15. Hz>, 'zz': <Quantity -6. Hz>}

        >>> b = PrincipalComponents(xx='2 Hz', yy='5 Hz', zz='-6 Hz')
        >>> print(b)
        {'xx': <Quantity 2. Hz>, 'yy': <Quantity 5. Hz>, 'zz': <Quantity -6. Hz>}

        >>> c = PrincipalComponents({'zz': '5 kHz'})
        >>> print(c)
        {'xx': <Quantity 0. Hz>, 'yy': <Quantity 0. Hz>, 'zz': <Quantity 5. kHz>}

    :return: An instance of :ref:`principal_api` class.
    """

    __slots__ = (
        '_dictionary_PAS',
    )

    def __init__(self, xx='0 Hz', yy='0 Hz', zz='0 Hz'):
        """Initialize an instance of PrincipalComponents class."""
        self._dictionary_PAS = {
            'xx': 0*u.Hz,
            'yy': 0*u.Hz,
            'zz': 0*u.Hz
        }

        self._dictionary_PAS = _initialize_class_attributes(
            self._dictionary_PAS, xx, yy, zz
        )
        l_ = self._dictionary_PAS['xx'].unit.physical_type
        m_ = self._dictionary_PAS['yy'].unit.physical_type
        r_ = self._dictionary_PAS['zz'].unit.physical_type
        if l_ != m_:
            raise Exception((
                "'xx' and 'yy' parameters must have the same physical type, found "
                f"{l_} and {m_}, respectively."
            ))
        if l_ != r_:
            raise Exception((
                "'xx' and 'zz' parameters must have the same physical type, found "
                f"{l_} and {r_}, respectively."
            ))
        if m_ != r_:
            raise Exception((
                "'yy' and 'zz' parameters must have the same physical type, found "
                f"{m_} and {r_}, respectively."
            ))

# private methods
    def __repr__(self):
        return f'<PrincipalComponents {str(self._dictionary_PAS)}>'

    def __str__(self):
        return str(self._dictionary_PAS)

# attributes
    @property
    def xx(self):
        r"""The x component of the tensor in the Principal axis system, :math:`xx`."""
        return self._dictionary_PAS['xx']

    @property
    def yy(self):
        r"""The y component of a tensor in the Principal axis system, :math:`yy`."""
        return self._dictionary_PAS['yy']

    @property
    def zz(self):
        r"""The z component of a tensor in the Principal axis system, :math:`zz`."""
        return self._dictionary_PAS['zz']

# methods
    def haeberlen_values(self):
        r"""
        The values in Haeberlen convension corresponding to the principal components.

        In Haeberlen convension the principal components of the second rank
        symmetric tensor is represented by three parameters, isotropic
        :math:`\delta`, anisotropic :math:`\zeta`, and asymmetry :math:`\eta`.

        The isotropic contribution is evaluated as

        :math:`\delta = (xx + yy + zz)/3`.

        The anisotropic, :math:`\zeta`, and the asymmetry, :math:`\eta`
        parameters are evaluated by first ordering the principle components
        such that

        :math:`|zz'| > |xx'| > |yy'|`

        where :math:`xx'=xx-\delta`, :math:`yy'=yy-\delta` and
        :math:`zz'=z-\delta`.

        Then, :math:`\zeta=zz'` and :math:`\eta=(yy'-xx')/\zeta`.

            >>> print(a.haeberlen_values())
            {'iso': <Quantity 3.33333333 Hz>, 'zeta': <Quantity -11.66666667 Hz>, 'eta': 0.6}

            >>> print(b.haeberlen_values())
            {'iso': <Quantity 0.33333333 Hz>, 'zeta': <Quantity -6.33333333 Hz>, 'eta': 0.4736842105263158}

            >>> print(c.haeberlen_values())
            {'iso': <Quantity 1666.66666667 Hz>, 'zeta': <Quantity 3.33333333 kHz>, 'eta': 0.0}

        :return: An instance of :ref:`haeberlen_api` class.
        """
        iso, zeta, eta = _get_haeberlen_values_from_PAS(self._dictionary_PAS)
        return HaeberlenConvension(iso=iso, zeta=zeta, eta=eta)

    def __dict_obj__(self):
        d_ = _get_json_object(self._dictionary_PAS)
        return d_

    def dumps(self):
        """
        Return the instance of the class as a json serialized string.

            >>> print(a.dumps())
            {
                "xx": "1.0 Hz",
                "yy": "15.0 Hz",
                "zz": "-6.0 Hz"
            }

        """
        d_ = _get_json_object(self._dictionary_PAS)
        return json.dumps(d_, indent=4, allow_nan=False, ensure_ascii=True)


class EulerAngles:
    r"""
    An instance of this class represents the three Euler angles, alpha, beta and gamma.

    The instance may be initialized using any of the following-

        >>> a = EulerAngles('3.12 rad', '3 deg')
        >>> print(a)
        {'alpha': <Quantity 3.12 rad>, 'beta': <Quantity 3. deg>, 'gamma': <Quantity 0. rad>}

        >>> b = EulerAngles(alpha='10 deg', beta='54.74 deg')
        >>> print(b)
        {'alpha': <Quantity 10. deg>, 'beta': <Quantity 54.74 deg>, 'gamma': <Quantity 0. rad>}

        >>> c = EulerAngles({'alpha':'0.1 rad','gamma':'45.2 deg'})
        >>> print(c)
        {'alpha': <Quantity 0.1 rad>, 'beta': <Quantity 0. rad>, 'gamma': <Quantity 45.2 deg>}

    :return: An instance of a :ref:`euler_angle_api` class.
    """

    __slots__ = (
        '_dictionary_euler_angle',
    )

    def __init__(self, alpha='0 rad', beta='0 rad', gamma='0 rad'):
        """Initialize an instance of EulerAngles class."""
        self._dictionary_euler_angle = {
            'alpha': 0*u.rad,
            'beta': 0*u.rad,
            'gamma': 0*u.rad
        }

        self._dictionary_euler_angle = _initialize_class_attributes(
            self._dictionary_euler_angle, alpha, beta, gamma
        )

        l_ = self._dictionary_euler_angle['alpha'].unit.physical_type
        m_ = self._dictionary_euler_angle['beta'].unit.physical_type
        r_ = self._dictionary_euler_angle['gamma'].unit.physical_type

        if l_ != m_:
            raise Exception((
                "'alpha' and 'beta' parameters must have the same physical type, found "
                f"{l_} and {m_}, respectively."
            ))
        if l_ != r_:
            raise Exception((
                "'alpha' and 'gamma' parameters must have the same physical type, found "
                f"{l_} and {r_}, respectively."
            ))
        if m_ != r_:
            raise Exception((
                "'beta' and 'gamma' parameters must have the same physical type, found "
                f"{m_} and {r_}, respectively."
            ))

# private methods
    def __repr__(self):
        return f'<EulerAngles {str(self._dictionary_euler_angle)}>'

    def __str__(self):
        return str(self._dictionary_euler_angle)

# attributes
    @property
    def alpha(self):
        """Return alpha Eular angle as a Quantity instance."""
        return self._dictionary_euler_angle['alpha']

    @property
    def beta(self):
        """Return beta Eular angle as a Quantity instance."""
        return self._dictionary_euler_angle['beta']

    @property
    def gamma(self):
        """Return gamma Eular angle as a Quantity instance."""
        return self._dictionary_euler_angle['gamma']

    # @property
    # def cartesian_coordinates(self):
    #     a_ = np.arccos(self.beta)
    #     return a

# methods
    def __dict_obj__(self):
        d_ = _get_json_object(self._dictionary_euler_angle)
        return d_

    def dumps(self):
        """
        Return the instance as a json serialized string.

            >>> print(a.dumps())
            {
                "alpha": "3.12 rad",
                "beta": "3.0 deg",
                "gamma": "0.0 rad"
            }

        """
        d_ = _get_json_object(self._dictionary_euler_angle)
        return json.dumps(d_, indent=4, allow_nan=False, ensure_ascii=True)
