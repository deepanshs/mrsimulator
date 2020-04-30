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
    Base isotopmer class representing an isolated spin-system with sites and couplings.

    Arguments:
        name: An optional string with the isotopomer name. The default is an empty
                string.
        description: An optional string describing the isotopomer. The default is
                an empty string.
        sites: A list of :ref:`site_api` objects or an equivalent python dict object
                representing a nuclear site. Default value is an empty list.
        abundance: The abundance of the isotopomer in unit of %. The default value is
                100. This attribute is useful when multiple isotopomers are present.
        transitions: A list of lists where each list is a transition, given as
                [initial state, final state]. For examples [[-0.5, 0.5]], [[1.5, 2.5]].
                You may also provide, multiple transitions, [[-0.5, 0.5], [-1.5, 0.5]].
                The resulting spectrum is a sum of resonances from individual
                transitions.
    """

    name: Optional[str] = ""
    description: Optional[str] = ""
    sites: List[Site] = []
    # couplings: list = [], # TODO: Deepansh what should this look like?
    abundance: float = Field(default=100, ge=0, le=100)

    property_unit_types: ClassVar = {"abundance": "dimensionless"}
    property_default_units: ClassVar = {"abundance": "pct"}
    property_units: Dict = {"abundance": "pct"}
    transitions: Optional[List] = None

    class Config:
        validate_assignment = True

    @classmethod
    def parse_dict_with_units(cls, py_dict):
        """
        Parse physical quantity from the attributes of an isotopomer object
        expressed as a python dictionary. The physical quantities are expressed
        as string with a number followed by a unit.

        Args:
            py_dict: Python dictionary representation of an isotopomer object with
                    attributes values as string with a physical quantity.

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

    def _Zeeman_energy_states(self):
        """
        Return the energy states as a Numpy array where axis 0 is the number of energy
        states and axis 1 is the spin quantum numbers. The spin quantum numbers are
        ordered based on the order of sites in the isotopomer.
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
    def Zeeman_energy_states(self):
        """
        Return a list of all Zeeman energy states, where the energy states are
        represented by a list of quantum numbers, [m_1, m_2,..m_n], where m_i is the
        quantum number associated with the ith site of the isotopomer.
        """
        states = self._Zeeman_energy_states()
        return [ZeemanState(len(self.sites), *item) for item in states]

    def _all_transitions(self):
        """
        Return transition as a Numpy array where axis 0 is the number of transitions,
        axis 1 is the initial energy state, and axis 2 is the final energy state.
        """
        energy_states = self._Zeeman_energy_states()
        s = energy_states.shape[0]
        lst = np.arange(s)
        indexes = np.asarray(np.meshgrid(lst, lst)).T.reshape(-1, 2)
        return energy_states[indexes]

    @property
    def all_transitions(self):
        """Returns a list of all possible spin transitions in as isotopomer."""
        transitions = self._all_transitions()
        return TransitionList(
            [
                Transition(initial=item[0].tolist(), final=item[1].tolist())
                for item in transitions
            ]
        )

    def to_freq_dict(self, B0):
        """
        Serialize the Isotopomer object to a JSON compliant python dictionary where the
        attribute values are numbers expressed in default units. The default unit
        for attributes with respective dimensionalities are:
        - frequency: `Hz`
        - angle: `rad`

        Args:
            B0: The macroscopic magnetic flux density in units of Tesla, T.

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
             'transitions': None}
        """
        temp_dict = self.dict()
        temp_dict["sites"] = [site.to_freq_dict(B0) for site in self.sites]
        temp_dict.pop("property_units")
        return temp_dict

    def get_isotopes(self, I=None):
        """
        An ordered list of isotopes from each site in an isotopomer corresponding to
        the given value of spin quantum number `I`. If `I` is None, a list of all
        isotopes is returned instead.

        Args:
            I: (optional) The spin quantum number. Valid input are multiples of 0.5.

        Returns:
            A list of isotopes.

        Example:
            >>> set(isotopomers.get_isotopes()) # doctest:+SKIP
            {'1H', '27Al', '13C'}
            >>> set(isotopomers.get_isotopes(I=0.5)) # doctest:+SKIP
            {'1H', '13C'}
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


def allowed_isotopes(I=None):
    """
    List of NMR active isotopes currently supported in ``mrsimulator``.

    Args:
        I: (optional) The spin quantum number. Valid input are multiples of 0.5.

    Returns:
        A list of all isotopes supported in ``mrsimulator`` with the given spin
        quantum number `I`. If the spin is unspecified or None, a list of all
        allowed isotopes is returned instead.
    Returns:
        Set
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
