# -*- coding: utf-8 -*-
"""Base Isotopomer class."""
from copy import deepcopy
from typing import ClassVar
from typing import Dict
from typing import List
from typing import Optional

import numpy as np
from pydantic import Field

from .abstract_list import TransitionList
from .isotope import ISOTOPE_DATA
from .parseable import Parseable
from .site import Site
from .state import ZeemanState
from .transition import Transition

__author__ = "Deepansh J. Srivastava"
__email__ = "deepansh2012@gmail.com"


class Isotopomer(Parseable):
    """
    Base isotopomer class representing an isolated spin-system containing multiple
    sites and couplings.

    Attributes:
        name: An optional string with the isotopomer name. The default is an empty
            string.
        description: An optional string describing the isotopomer. The default is
            an empty string.
        sites: A list of :ref:`site_api` objects or a list of equivalent python
            dictionary objects representing an NMR site. The default value is an empty
            list.
        abundance: The abundance of isotopomer in the units of %. The default value
            is 100. This attribute is useful when multiple isotopomers are present.
        transition_pathways: A list of lists, where the inner list represents a
            transition pathway, and the outer list is the number of transition
            pathways. Each transition pathways is a list of Transition objects.

            An example of transition pathways follows

            .. doctest::

                >>> isotopomer_1.transition_pathways = [
                ...   [{'initial': -0.5, 'final': 0.5}, {'initial': 1.5, 'final': 2.5}]
                ... ]

            The resulting spectrum is a sum of resonances arising from individual
            transition pathways.

            .. note::
                The `transition_pathways` is a transient attribute that overrides the
                transition pathway query from the NMR method object. Only use this
                attribute if you know what you are doing. This attribute may provide
                a significant improvement in the performance in iterative algorithms,
                such as the least-squares algorithm, where a one-time query is
                sufficient. Remember, repeated queries for the transition pathways will
                add significant overhead to the computation.
    """

    name: Optional[str] = ""
    description: Optional[str] = ""
    sites: List[Site] = []
    # couplings: list = [], # TODO: Deepansh what should this look like?
    abundance: float = Field(default=100, ge=0, le=100)

    property_unit_types: ClassVar = {"abundance": "dimensionless"}
    property_default_units: ClassVar = {"abundance": "pct"}
    property_units: Dict = {"abundance": "pct"}
    transition_pathways: Optional[List] = None

    class Config:
        validate_assignment = True

    def _Zeeman_energy_states(self) -> np.ndarray:
        """
        Return the energy states as a Numpy array where the axis 0 is the number of
        energy states, and axis 1 is the spin quantum numbers. The spin quantum numbers
        are ordered based on the order of the sites within the isotopomer.
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

    @property
    def Zeeman_energy_states(self) -> list:
        r"""
        Return a list of all Zeeman energy states of the isotopomer spin-system,
        where the energy states are represented by a list of quantum numbers,

        .. math::
            |\text{state}⟩ = [m_1, m_2,.. m_n],

        where :math:`m_i` is the quantum number associated with the :math:`i^\text{th}`
        site within the isotopomer.

        Example:
            >>> isotopomer_1H_13C.get_isotopes() # two site (spin-1/2) isotopomer
            ['13C', '1H']
            >>> isotopomer_1H_13C.Zeeman_energy_states  # four energy level system.
            [|-0.5, -0.5⟩, |-0.5, 0.5⟩, |0.5, -0.5⟩, |0.5, 0.5⟩]

        Return: A list of ZeemanState objects.
        """
        states = self._Zeeman_energy_states()
        return [ZeemanState(len(self.sites), *item) for item in states]

    def _all_transitions(self) -> np.ndarray:
        """
        Return all transitions of an isotopomer as a Numpy array of shape (M, 2, N),
        where M is the number of transitions, and N is the number of sites in the
        isotopomer. The second axis is of length 2, where the entries at T[:, 0, :]
        are the initial energy states, and the entries at T[:, 1, :], the corresponding
        final energy states of the spin transitions.
        """
        energy_states = self._Zeeman_energy_states()
        s = energy_states.shape[0]
        lst = np.arange(s)
        indexes = np.asarray(np.meshgrid(lst, lst)).T.reshape(-1, 2)
        return energy_states[indexes]

    @property
    def all_transitions(self) -> TransitionList:
        """Returns a list of all possible spin transitions in the given isotopomer.

        Example:
            >>> isotopomer_1H_13C.get_isotopes()  # two site (spin-1/2) isotopomer
            ['13C', '1H']
            >>> isotopomer_1H_13C.all_transitions  # 16 two energy level transitions
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
        Parse the physical quantity from the attributes of an isotopomer object
        when represented as a python dictionary. The physical quantities are
        expressed as a string with a number followed by a unit.

        Args:
            dict py_dict: A python dictionary representation of an isotopomer object
                where attributes values are given as a string with a physical quantity.

        Example:
            >>> isotopomer_dict = {
            ...     "sites": [{
            ...         "isotope":"13C",
            ...         "isotropic_chemical_shift": "20 ppm",
            ...         "shielding_symmetric": {
            ...             "zeta": "10 ppm",
            ...             "eta": 0.5
            ...         }
            ...     }]
            ... }
            >>> isotopomer_1 = Isotopomer.parse_dict_with_units(isotopomer_dict)
        """
        py_dict_copy = deepcopy(py_dict)
        if "sites" in py_dict_copy:
            py_dict_copy["sites"] = [
                Site.parse_dict_with_units(s) for s in py_dict_copy["sites"]
            ]

        return super().parse_dict_with_units(py_dict_copy)

    def to_freq_dict(self, B0: float) -> dict:
        """
        Serialize the Isotopomer object to a JSON compliant python dictionary object
        where the attribute values are numbers expressed in default units. The default
        unit for the attributes with respective dimensionalities are:

        - frequency: `Hz`
        - angle: `rad`

        Args:
            float B0: The macroscopic magnetic flux density in units of Tesla, T.

        Return: A python dict

        Example:
            >>> pprint(isotopomer_1.to_freq_dict(B0=9.4))
            {'abundance': 100,
             'description': '',
             'name': '',
             'sites': [{'isotope': '13C',
                        'isotropic_chemical_shift': -2013.1791999999998,
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

    def get_isotopes(self, I=None) -> list:
        """
        An ordered list of isotopes from sites within the isotopomer corresponding to
        the given value of spin quantum number `I`. If `I` is None, a list of all
        isotopes is returned instead.

        Args:
            float I: An optional spin quantum number. The valid inputs are the
                multiples of 0.5.

        Returns:
            A list of isotopes.

        Example:
            >>> isotopomers.get_isotopes()
            ['13C', '1H', '27Al']
            >>> isotopomers.get_isotopes(I=0.5)
            ['13C', '1H']
            >>> isotopomers.get_isotopes(I=1.5)
            []
            >>> isotopomers.get_isotopes(I=2.5)
            ['27Al']
        """
        isotope_list = allowed_isotopes(I)
        return [
            site.isotope.symbol
            for site in self.sites
            if site.isotope.symbol in isotope_list
        ]


def allowed_isotopes(I=None) -> list:
    """
    List of NMR active isotopes currently supported in ``mrsimulator``.

    Args:
        I: An optional float with the spin quantum number. The valid inputs are the
            multiples of 0.5.

    Returns:
        A list of all isotopes supported in ``mrsimulator`` with the given spin
        quantum number `I`. If the spin is unspecified or None, a list of all
        allowed isotopes is returned instead.
    """
    if I is None:
        return list(ISOTOPE_DATA.keys())
    return list(
        {
            isotope
            for isotope, data in ISOTOPE_DATA.items()
            if data["spin"] == int(2 * I)
        }
    )
