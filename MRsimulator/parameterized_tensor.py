
r"""
The ParameterizedTensor class is a collection of parmeterized nuclear
shielding and electric quadrupole tensor with parameters that are commonly
encountered in the Nuclear Magnetic Resonance (NMR).
"""


from astropy import units as u
from copy import deepcopy
import json
from .core import (
    HaeberlenConvension,
    PrincipalComponents,
    EulerAngles,
    _initialize_class_attributes,
    _get_json_object
)
from warnings import warn

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.90@osu.edu"


def _sub_class(func, value):
    if isinstance(value, list) or isinstance(value, tuple):
        value = func(*value)
    else:
        value = func(value)
    return value

def _get_kwargs(args, kwargs, class_name):
    # print(args)
    # print(kwargs)
    if len(args) != 0 and kwargs != {}:
        raise Exception((
                f"Cannot provide both argument, {list(args)} and "
                f"keyword arguments {list(kwargs)} at the same time "
                f"for {class_name}."
                )
            )
    if len(args) == 1:
        if isinstance(args[0], dict):
            kwargs = args[0]
            return kwargs
        if isinstance(args[0], list) or isinstance(args[0], tuple):
            if class_name == 'ElectricQuadrupoleTensor':
                if len(args[0]) > 3:
                    raise Exception((
                        f"Failed to create a {class_name} instance. Cannot "
                        f"interpret the arguments {args[0]}."))
                kwargs = {}
                list__ = ['Cq', 'eta', 'euler_angles']
                for i in range(len(args[0])):
                    kwargs[list__[i]] = args[0][i]
                warn__ = (
                    f"When using a list or tuple of values, the values "
                    f"will be assigned as 'Cq', 'eta' and 'euler_angles' "
                    f"respectively."
                )
                warn(str(warn__))
                return kwargs
            
            if class_name == 'NuclearShieldingTensor':
                if len(args[0]) > 2:
                    raise Exception((
                        f"Failed to create a {class_name} instance. Cannot "
                        f"interpret the arguments {args[0]}."))
                kwargs = {}
                list__ = ['harberlen', 'euler_angles']
                for i in range(len(args[0])):
                    kwargs[list__[i]] = args[0][i]
                warn__ = (
                    f"When using a list or tuple of values, the ordered "
                    f"values will be assigned as 'haeberlen' and 'euler_angles' "
                    f"respectively."
                )
                warn(str(warn__))
                return kwargs

        else:
            raise Exception((
                "Only a dictionary of valid keywords is allowed as an "
                f"argument of {class_name}, found {type(args[0])} with "
                f"value '{args[0]}'."
                )
            )
    elif len(args) > 1:
        raise Exception((
            f"When passing arguments to {class_name}, only one argument of "
            f"type 'dict' is allowed, {len(args)} arguments are provided."
        )
    )
    return kwargs



class ParameterizedTensor:
    """Abstract class"""

    __slots__ = (
        '_euler_angles',
        '_principal_components',
        '_haeberlen'
    )
    def __init__(self, 
                 euler_angles=EulerAngles(),
                 principal_components=PrincipalComponents(),
                 haeberlen=HaeberlenConvension()):

        self._euler_angles = euler_angles
        self._principal_components = principal_components
        self._haeberlen = haeberlen

# attributes
    @property
    def euler_angles(self):
        r"""
        Return an instance of :ref:`euler_angle_api`.

        The Euler angles are the set of three angles, :math:`\alpha`,
        :math:`\beta` and :math:`\gamma` that tranforms the principal
        components of the second rank symmetric tensor to a common
        crystal/molecular frame.
        """
        return self._euler_angles

    @property
    def principal_components(self):
        r"""
        Return an instance of :ref:`principal_api`.

        The components, :math:`T_{xx}`, :math:`T_{yy}` and :math:`T_{zz}` of the
        second rank symmetric tensor, :math:`T`, in the principal axis system.
        """
        return self._principal_components

    @property
    def haeberlen(self):
        r"""
        Return an instance of :ref:`haeberlen_api` class.

        The attributes of the HaeberlenConvension instance are the
        nuclear shielding tensor principal conponents represented in
        Haeberlen convension.
        """
        return self._haeberlen

# methods
    def __dict_obj__(self):
        d_ = {key[1:]: self.__getattribute__(key) for key in self.__class__.__list__}
        return d_

    def dumps(self):
        """
        Return the instance as a json serialized string.
        
            >>> from mrsimulator.parameterized_tensor import NuclearShieldingTensor
            >>> n1 = NuclearShieldingTensor(haeberlen={'iso': '1210 Hz', 'zeta': '10 kHz', 'eta': 0.1}, euler_angles={'alpha':'5 deg','gamma':'54.74 deg'})
            >>> print(n1.dumps())
            {
                "haeberlen": {
                    "iso": "1210.0 Hz",
                    "zeta": "10.0 kHz",
                    "eta": 0.1
                },
                "euler_angles": {
                    "alpha": "5.0 deg",
                    "beta": "0.0 rad",
                    "gamma": "54.74 deg"
                }
            }
        """
        # d_ = {key[1:]: self.__getattribute__(key) for key in self.__class__.__list__}
        d_ = _get_json_object(self.__dict_obj__())
        return json.dumps(d_, indent=4, allow_nan=False, ensure_ascii=True)


class ElectricQuadrupoleTensor(ParameterizedTensor):
    r"""
    The class describes a parameterized second rank symmetric electric
    quadrupole tensor. There are several ways to initialize an instance of
    the ElectricQuadrupoleTensor class.

    * **Initialize with** :math:`C_q`, :math:`\eta` **and** :math:`\Theta`.
    
    Here, :math:`C_q` is the quadrupolar coupling constant, :math:`\eta` is the
    quadrupolar asymmetry parameter, and :math:`\Theta = [\alpha, \beta, \gamma]`
    are the three euler angles. 
    The three parameters may be specified through keyword arguments as
    
    .. doctest::

        >>> from mrsimulator.parameterized_tensor import ElectricQuadrupoleTensor

        >>> q_keywords = ElectricQuadrupoleTensor(
        ...         Cq='10 MHz',
        ...         eta=1,
        ...         euler_angles={'alpha':'5 deg','gamma':'54.74 deg'}
        ...     )
        >>> print(q_keywords.dumps())
        {
            "Cq": "10.0 MHz",
            "eta": 1.0,
            "euler_angles": {
                "alpha": "5.0 deg",
                "beta": "0.0 rad",
                "gamma": "54.74 deg"
            }
        }

    as a dictionary of keywords,

    .. doctest ::

        >>> q_dict = ElectricQuadrupoleTensor(
        ...         {   
        ...             'Cq': '5.1 MHz',
        ...             'eta': 0.5,
        ...             'euler_angles':('5 deg','54.74 deg')
        ...         }
        ...     )
        >>> print(q_dict.dumps())
        {
            "Cq": "5.1 MHz",
            "eta": 0.5,
            "euler_angles": {
                "alpha": "5.0 deg",
                "beta": "54.74 deg",
                "gamma": "0.0 rad"
            }
        }

    as a list/tuple of values,

    .. doctest ::

        >>> q_list = ElectricQuadrupoleTensor( ('5.1 MHz', 0.5, ('5 deg','54.74 deg')) )
        >>> print(q_list.dumps())
        {
            "Cq": "5.1 MHz",
            "eta": 0.5,
            "euler_angles": {
                "alpha": "5.0 deg",
                "beta": "54.74 deg",
                "gamma": "0.0 rad"
            }
        }

    .. note::
        The value of ``euler_angles`` keyword is forwarded to the EulerAngles class.
        For more information refer to the :ref:`euler_angle_api` class. 

    The default value of the ElectricQuadrupoleTensor instance is

        >>> q_default = ElectricQuadrupoleTensor()
        >>> print(q_default.dumps())
        {
            "Cq": "0.0 Hz",
            "eta": 0.0,
            "euler_angles": {
                "alpha": "0.0 rad",
                "beta": "0.0 rad",
                "gamma": "0.0 rad"
            }
        }

    :return: An instance of the :ref:`quad_tensor_api` class.
    """
    __slots__ = (
        '_dictionary_quad',
        '_Cq',
        '_eta'
    )

    __list__ = ['_Cq', '_eta', '_euler_angles']

    def __new__(cls, *args, **kwargs):
        if len(args) == 1:
            if isinstance(args[0], cls):
                return args[0]
        instance = super(ElectricQuadrupoleTensor, cls).__new__(cls)
        return instance

    def __init__(self, *args, **kwargs):
        """Initialize an instance of ElectricQuadrupoleTensor class."""

        _dictionary_quad = {
            'Cq': 0*u.Hz,
            'eta': 0.0,
            'euler_angles': EulerAngles(),
            'principal_components': PrincipalComponents()
        }

        kwargs = _get_kwargs(args, kwargs, self.__class__.__name__)
        
        keys = kwargs.keys()
        if 'euler_angles' in keys:
            kwargs['euler_angles'] = _sub_class(
                EulerAngles, kwargs['euler_angles']
            )

        _dictionary_quad = _initialize_class_attributes(
            _dictionary_quad, **kwargs
        )

        super().__init__(
            _dictionary_quad['euler_angles'],
            _dictionary_quad['principal_components']
        )
        self._Cq = _dictionary_quad['Cq']
        self._eta = _dictionary_quad['eta']

# private methods
    def __repr__(self):
        return f'<ElectricQuadrupoleTensor {str(self.__dict_obj__())}>'

    def __str__(self):
        return str(self.__dict_obj__())

# attributes
    @property
    def Cq(self):
        r"""The electric quadrupole coupling constant."""
        return self._Cq #_dictionary_quad["Cq"]

    @property
    def eta(self):
        r"""The quadrupole asymmetry parameter."""
        return self._eta #_dictionary_quad["eta"]


class NuclearShieldingTensor(ParameterizedTensor):
    r"""
    The class describes a parameterized second rank symmetric nuclear shielding
    tensor. There are several ways to initialize an instance of
    the NuclearShieldingTensor class.

    * **Initialize with haeberlen values and euler angle**.

    The parameters may be specified through keyword arguments as

    .. doctest::

        >>> from mrsimulator.parameterized_tensor import NuclearShieldingTensor

        >>> n_keywords = NuclearShieldingTensor(
        ...         haeberlen={'iso': '1210 Hz', 'zeta': '10 kHz', 'eta': 0.1},
        ...         euler_angles={'alpha':'5 deg','gamma':'54.74 deg'}
        ...     )
        >>> print(n_keywords.dumps())
        {
            "haeberlen": {
                "iso": "1210.0 Hz",
                "zeta": "10.0 kHz",
                "eta": 0.1
            },
            "euler_angles": {
                "alpha": "5.0 deg",
                "beta": "0.0 rad",
                "gamma": "54.74 deg"
            }
        }

    or as a dictionary of keywords,

    .. doctest::
    
        >>> n_dict = NuclearShieldingTensor(
        ...         {
        ...             'haeberlen': {'iso': '1210 Hz', 'zeta': '10 kHz', 'eta': 0.1},
        ...             'euler_angles': {'alpha':'5 deg','gamma':'54.74 deg'}
        ...         }
        ...     )
        >>> print(n_dict.dumps())
        {
            "haeberlen": {
                "iso": "1210.0 Hz",
                "zeta": "10.0 kHz",
                "eta": 0.1
            },
            "euler_angles": {
                "alpha": "5.0 deg",
                "beta": "0.0 rad",
                "gamma": "54.74 deg"
            }
        }


    * **Initialize with principal components and euler angle**.

    As keywords..

    .. doctest::

        >>> n_keywords = NuclearShieldingTensor(
        ...         principal_components=['-5000 Hz', '-5000 Hz', '10 kHz'],
        ...         euler_angles=('5 deg', '54.74 deg')
        ...     )
        >>> print(n_keywords.dumps())
        {
            "haeberlen": {
                "iso": "0.0 Hz",
                "zeta": "10.0 kHz",
                "eta": 0.0
            },
            "euler_angles": {
                "alpha": "5.0 deg",
                "beta": "54.74 deg",
                "gamma": "0.0 rad"
            }
        }

    or as a dictionary of keywords,

    .. doctest::

        >>> n_dict = NuclearShieldingTensor(
        ...         {
        ...             'principal_components': {'xx': '-5000 Hz', 'yy': '-1000 Hz', 'zz': '10 kHz'},
        ...             'euler_angles': {'alpha':'5 deg','gamma':'54.74 deg'}
        ...         }
        ...     )
        >>> print(n_dict.dumps())
        {
            "haeberlen": {
                "iso": "1333.3333333333333 Hz",
                "zeta": "8.666666666666666 kHz",
                "eta": 0.46153846153846156
            },
            "euler_angles": {
                "alpha": "5.0 deg",
                "beta": "0.0 rad",
                "gamma": "54.74 deg"
            }
        }

    .. note::
        The value of ``euler_angles``, ``haeberlen`` and ``principal_components``
        keywords are forwarded to :ref:`euler_angle_api`, :ref:`haeberlen_api`
        and :ref:`principal_api` classes respectively.
        For more information refer to the documentation of these class. 

    The default value of the NuclearShieldingTensor instance is

        >>> n_default = NuclearShieldingTensor()
        >>> print(n_default.dumps())
        {
            "haeberlen": {
                "iso": "0.0 Hz",
                "zeta": "0.0 Hz",
                "eta": 0.0
            },
            "euler_angles": {
                "alpha": "0.0 rad",
                "beta": "0.0 rad",
                "gamma": "0.0 rad"
            }
        }

    :return: An instance of the :ref:`shielding_tensor_api` class.
    """

    __slots__ = (
        '_dictionary_nuclear_shielding',
        '_haeberlen',
    )

    __list__ = ['_haeberlen', '_euler_angles']

    def __new__(cls, *args, **kwargs):
        if len(args) == 1:
            if isinstance(args[0], cls):
                return args[0]
        instance = super(NuclearShieldingTensor, cls).__new__(cls)
        return instance

    def __init__(self, *args, **kwargs):
        """Instantiate NuclearShieldingTensor class."""

        _dictionary_nuclear_shielding = {
            'haeberlen': HaeberlenConvension(),
            'principal_components': PrincipalComponents(),
            'euler_angles': EulerAngles()
        }

        kwargs = _get_kwargs(args, kwargs, self.__class__.__name__)

        keys = list(kwargs.keys())
        if 'euler_angles' in keys:
            kwargs['euler_angles'] = _sub_class(EulerAngles, kwargs['euler_angles'])

        if 'haeberlen' in keys:
            kwargs['haeberlen'] = _sub_class(HaeberlenConvension, kwargs['haeberlen'])
            kwargs['principal_components'] = kwargs['haeberlen'].principal_values()

        if 'principal_components' in keys:
            kwargs['principal_components'] = _sub_class(
                PrincipalComponents, kwargs['principal_components']
            )
            kwargs['haeberlen'] = kwargs['principal_components'].haeberlen_values()

        _dictionary_nuclear_shielding = _initialize_class_attributes(
            _dictionary_nuclear_shielding, **kwargs
        )

        super().__init__(
            _dictionary_nuclear_shielding['euler_angles'],
            _dictionary_nuclear_shielding['principal_components'],
            _dictionary_nuclear_shielding['haeberlen']
        )

# private methods
    def __repr__(self):
        return f'<NuclearShieldingTensor {str(self.__dict_obj__())}>'

    def __str__(self):
        return f'{str(self.__dict_obj__())}'

    def __dict_obj__(self):
        d_ = {key[1:]: self.__getattribute__(key) for key in self.__class__.__list__}
        return d_



