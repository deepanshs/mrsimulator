# -*- coding: utf-8 -*-
"""Base SpinSystem class."""
from copy import deepcopy
from typing import ClassVar
from typing import Dict
from typing import List
from typing import Optional

import numpy as np
from mrsimulator.transition import Transition
from mrsimulator.transition.transition_list import TransitionList
from mrsimulator.utils.parseable import Parseable
from pydantic import Field

from .isotope import ISOTOPE_DATA
from .site import Site
from .zeeman_state import ZeemanState

__author__ = "Deepansh J. Srivastava"
__email__ = "deepansh2012@gmail.com"


class SpinSystem(Parseable):
    """
    Base class representing an isolated spin system containing multiple sites and
    couplings amongst them.

    .. rubric:: Attribute Documentation

    Attributes
    ----------

    sites: A list of :ref:`site_api` objects or equivalent dict objects (optional).
        The value is a list of sites within the spin system, where each site represents
        a single-site nuclear spin interaction tensor parameters. The default value is
        an empty list.

        Example
        -------

        >>> sys1 = SpinSystem()
        >>> sys1.sites = [Site(isotope='17O'), Site(isotope='1H')]
        >>> # or equivalently
        >>> sys1.sites = [{'isotope': '17O'}, {'isotope': '1H'}]

    abundance: float (optional).
        The abundance of the spin system in units of %. The default value is 100. The
        value of this attribute is useful when multiple spin systems are present.

        Example
        -------

        >>> sys1.abundance = 10

    name: str (optional).
        The value is the name or id of the spin system. The default value is None.

        Example
        -------

        >>> sys1.name = '1H-17O-0'
        >>> sys1.name
        '1H-17O-0'

    label: str (optional).
        The value is a label for the spin system. The default value is None.

        Example
        -------

        >>> sys1.label = 'Heteronuclear spin system'
        >>> sys1.label
        'Heteronuclear spin system'

    description: str (optional).
        The value is a description of the spin system. The default value is None.

        Example
        -------

        >>> sys1.description = 'A test for the spin system'
        >>> sys1.description
        'A test for the spin system'

    transition_pathways: list (optional).
        The value is a list of lists, where the inner list represents a transition
        pathway, and the outer list is the number of transition pathways. Each
        transition pathways is a list of Transition objects. The resulting spectrum is
        a sum of the resonances arising from individual transition pathways. The
        default value is None.

        Example
        -------
        >>> sys1.transition_pathways = [
        ...     [
        ...         {'initial': [-2.5, 0.5], 'final': [2.5, 0.5]},
        ...         {'initial': [0.5, 0.5], 'final': [-0.5, 0.5]}
        ...     ]
        ... ]

        .. note::
            From any given spin system, the list of relevant transition pathways is
            determined by the applied NMR method. For example, consider a single site
            I=3/2 spin system. For this system, a Bloch decay spectrum method will
            select three transition pathways, one corresponding to the central and two
            to the satellite transitions. On the other hand, a Bloch decay central
            transition selective method will only select one transition pathway,
            corresponding to the central transition.

            Since the spin system is independent of the NMR method, the value of this
            attribute is, therefore, transient. You may use this attribute to override
            the default transition pathway query selection criterion of the NMR method
            objects.

            **Only use this attribute if you know what you are doing.**

            At times, this attribute may provide a significant improvement in the
            performance, especially in iterative algorithms, such as the least-squares
            algorithm, where a one-time transition pathway query is sufficient.
            Repeated queries for the transition pathways will add significant overhead
            to the computation.

            .. seealso::
                `Fitting example
                <./../auto_examples/Fitting/plot_2_mrsimFitExample_O17.html>`_
    """

    name: str = None
    label: str = None
    description: str = None
    sites: List[Site] = []
    # couplings: list = [], # TODO: Deepansh what should this look like?
    abundance: float = Field(default=100, ge=0, le=100)
    transition_pathways: Optional[List] = None

    property_unit_types: ClassVar = {"abundance": "dimensionless"}
    property_default_units: ClassVar = {"abundance": "pct"}
    property_units: Dict = {"abundance": "pct"}

    class Config:
        validate_assignment = True

    def get_isotopes(self, spin_I=None) -> list:
        """
        An ordered list of isotopes from sites within the spin system corresponding to
        the given value of spin quantum number `I`. If `I` is None, a list of all
        isotopes is returned instead.

        Args:
            float spin_I: An optional spin quantum number. The valid inputs are the
                multiples of 0.5.

        Returns:
            A list of isotopes.

        Example
        -------

        >>> spin_systems.get_isotopes() # three spin systems
        ['13C', '1H', '27Al']
        >>> spin_systems.get_isotopes(spin_I=0.5) # isotopes with I=0.5
        ['13C', '1H']
        >>> spin_systems.get_isotopes(spin_I=1.5) # isotopes with I=1.5
        []
        >>> spin_systems.get_isotopes(spin_I=2.5) # isotopes with I=2.5
        ['27Al']
        """
        isotope_list = allowed_isotopes(spin_I)
        return [
            site.isotope.symbol
            for site in self.sites
            if site.isotope.symbol in isotope_list
        ]

    def _zeeman_energy_states(self) -> np.ndarray:
        """
        Return the energy states as a Numpy array where the axis 0 is the number of
        energy states, and axis 1 is the spin quantum numbers. The spin quantum numbers
        are ordered based on the order of the sites within the spin system.
        """
        two_I_p_one = [int(2 * site.isotope.spin + 1) for site in self.sites]
        spin_quantum_numbers = [
            np.arange(2 * site.isotope.spin + 1) - site.isotope.spin
            for site in self.sites
        ]
        size = len(spin_quantum_numbers)

        lst = []
        for j in range(size):
            k = 1
            for i in range(size):
                if i == j:
                    k = np.kron(k, spin_quantum_numbers[i])
                else:
                    k = np.kron(k, np.ones(two_I_p_one[i]))
            lst.append(k)
        return np.asarray(lst).T

    def zeeman_energy_states(self) -> list:
        r"""
        Return a list of all Zeeman energy states of the spin system,
        where the energy states are represented by a list of quantum numbers,

        .. math::
            |\Psi⟩ = [m_1, m_2,.. m_n],

        where :math:`m_i` is the quantum number associated with the :math:`i^\text{th}`
        site within the spin system, and :math:`\Psi` is the energy state.

        Example
        -------

        >>> spin_system_1H_13C.get_isotopes() # two site (spin-1/2) spin systems
        ['13C', '1H']
        >>> spin_system_1H_13C.zeeman_energy_states()  # four energy level system.
        [|-0.5, -0.5⟩, |-0.5, 0.5⟩, |0.5, -0.5⟩, |0.5, 0.5⟩]

        Returns
            A list of :ref:`zeeman_api` objects.
        """
        states = self._zeeman_energy_states()
        return [ZeemanState(len(self.sites), *item) for item in states]

    def _all_transitions(self) -> np.ndarray:
        """
        Return all transitions from a spin system as a Numpy array of shape (M, 2, N),
        where M is the number of transitions, and N is the number of sites in the
        spin system. The second axis is of length 2, where the entries at T[:, 0, :]
        are the initial energy states, and the entries at T[:, 1, :], the corresponding
        final energy states of the spin transitions.
        """
        energy_states = self._zeeman_energy_states()
        s = energy_states.shape[0]
        lst = np.arange(s)
        indexes = np.asarray(np.meshgrid(lst, lst)).T.reshape(-1, 2)
        return energy_states[indexes]

    def all_transitions(self) -> TransitionList:
        """Returns a list of all possible spin transitions in the given spin system.

        Example
        -------

        >>> spin_system_1H_13C.get_isotopes()  # two site (spin-1/2) spin system
        ['13C', '1H']
        >>> spin_system_1H_13C.all_transitions()  # 16 two energy level transitions
        [|-0.5, -0.5⟩⟨-0.5, -0.5|,
        |-0.5, 0.5⟩⟨-0.5, -0.5|,
        |0.5, -0.5⟩⟨-0.5, -0.5|,
        |0.5, 0.5⟩⟨-0.5, -0.5|,
        |-0.5, -0.5⟩⟨-0.5, 0.5|,
        |-0.5, 0.5⟩⟨-0.5, 0.5|,
        |0.5, -0.5⟩⟨-0.5, 0.5|,
        |0.5, 0.5⟩⟨-0.5, 0.5|,
        |-0.5, -0.5⟩⟨0.5, -0.5|,
        |-0.5, 0.5⟩⟨0.5, -0.5|,
        |0.5, -0.5⟩⟨0.5, -0.5|,
        |0.5, 0.5⟩⟨0.5, -0.5|,
        |-0.5, -0.5⟩⟨0.5, 0.5|,
        |-0.5, 0.5⟩⟨0.5, 0.5|,
        |0.5, -0.5⟩⟨0.5, 0.5|,
        |0.5, 0.5⟩⟨0.5, 0.5|]
        """
        transitions = self._all_transitions()
        return TransitionList(
            [
                Transition(initial=item[0].tolist(), final=item[1].tolist())
                for item in transitions
            ]
        )

    @classmethod
    def parse_dict_with_units(cls, py_dict: dict) -> dict:
        """
        Parse the physical quantity from a dictionary representation of the SpinSystem
        object, where the physical quantity is expressed as a string with a number and
        a unit.

        Args:
            dict py_dict: A required python dict object.

        Returns:
            :ref:`spin_sys_api` object.

        Example
        -------

        >>> spin_system_dict = {
        ...     "sites": [{
        ...         "isotope":"13C",
        ...         "isotropic_chemical_shift": "20 ppm",
        ...         "shielding_symmetric": {
        ...             "zeta": "10 ppm",
        ...             "eta": 0.5
        ...         }
        ...     }]
        ... }
        >>> spin_system_1 = SpinSystem.parse_dict_with_units(spin_system_dict)
        """
        py_dict_copy = deepcopy(py_dict)
        if "sites" in py_dict_copy:
            py_dict_copy["sites"] = [
                Site.parse_dict_with_units(s) for s in py_dict_copy["sites"]
            ]

        return super().parse_dict_with_units(py_dict_copy)

    def to_freq_dict(self, B0: float) -> dict:
        """
        Serialize the SpinSystem object to a JSON compliant python dictionary object,
        where the attribute value is a numbers expressed in the attribute's default
        unit. The default unit for the attributes with respective dimensionalities are:

        - frequency: `Hz`
        - angle: `rad`

        Args:
            float B0: A required macroscopic magnetic flux density in units of T.

        Return:
            A python dict

        Example
        -------

        >>> pprint(spin_system_1.to_freq_dict(B0=9.4))
        {'abundance': 100,
         'description': None,
         'label': None,
         'name': None,
         'sites': [{'description': None,
                    'isotope': '13C',
                    'isotropic_chemical_shift': -2013.1791999999998,
                    'label': None,
                    'name': None,
                    'quadrupolar': None,
                    'shielding_antisymmetric': None,
                    'shielding_symmetric': {'alpha': None,
                                            'beta': None,
                                            'eta': 0.5,
                                            'gamma': None,
                                            'zeta': -1006.5895999999999}}],
         'transition_pathways': None}
        """
        temp_dict = self.dict()
        temp_dict["sites"] = [site.to_freq_dict(B0) for site in self.sites]
        temp_dict.pop("property_units")
        return temp_dict


def allowed_isotopes(spin_I=None) -> list:
    """
    List of NMR active isotopes currently supported in ``mrsimulator``.

    Args:
        spin_I: An optional float with the spin quantum number. The valid inputs are
            the multiples of 0.5.

    Returns:
        A list of all isotopes supported in ``mrsimulator`` with the given spin
        quantum number `I`. If the spin is unspecified or None, a list of all
        allowed isotopes is returned instead.
    """
    if spin_I is None:
        return list(ISOTOPE_DATA.keys())
    return list(
        {
            isotope
            for isotope, data in ISOTOPE_DATA.items()
            if data["spin"] == int(2 * spin_I)
        }
    )
