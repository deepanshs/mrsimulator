# -*- coding: utf-8 -*-
"""Base SpinSystem class."""
from copy import deepcopy
from typing import ClassVar
from typing import Dict
from typing import List
from typing import Union

import numpy as np
from mrsimulator.base_model import get_zeeman_states
from mrsimulator.transition import Transition
from mrsimulator.transition.pathway import TransitionList
from mrsimulator.transition.pathway import TransitionPathway
from mrsimulator.utils.parseable import Parseable
from pydantic import Field
from pydantic import validator

from .coupling import Coupling
from .isotope import ISOTOPE_DATA
from .site import Site
from .zeemanstate import ZeemanState

__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"


class SpinSystem(Parseable):
    """
    Base class representing an isolated spin system containing multiple sites and
    couplings amongst them.

    .. rubric:: Attribute Documentation

    Attributes
    ----------

    sites: list of :ref:`site_api` or equivalent dict objects (optional)
        A list of :ref:`site_api` or equivalent dict objects within the spin system.
        Each site object represents single-site nuclear spin interaction (nuclear
        shielding and EFG) tensor parameters. The default value is an empty list.

        Example
        -------

        >>> sys1 = SpinSystem()
        >>> sys1.sites = [Site(isotope="17O"), Site(isotope="1H")]
        >>> # or equivalently
        >>> sys1.sites = [{"isotope": "17O"}, {"isotope": "1H"}]

    couplings: list of :ref:`coupling_api` or equivalent dict objects (optional)
        A list of :ref:`coupling_api` or equivalent dict objects within the spin
        system. Each coupling object represents two-site spin interaction (J-coupling
        and Dipolar) tensor parameters. The default value is an empty list.

        Example
        -------

        >>> sys1 = SpinSystem()
        >>> sys1.couplings = [
        ...     Coupling(site_index=[0, 1], isotropic_j=10.1),
        ...     Coupling(site_index=[2, 1], dipolar={"D": 1500})
        ... ]
        >>> # or equivalently
        >>> sys1.couplings = [
        ...     {"site_index": [0, 1], "isotropic_j": 10.1},
        ...     {"site_index": [2, 1], "dipolar": {"D": 1500}}
        ... ]

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

        >>> sys1.name = "1H-17O-0"
        >>> print(sys1.name)
        1H-17O-0

    label: str (optional).
        The value is a label for the spin system. The default value is None.

        Example
        -------

        >>> sys1.label = "Heteronuclear spin system"
        >>> print(sys1.label)
        Heteronuclear spin system

    description: str (optional).
        The value is a description of the spin system. The default value is None.

        Example
        -------

        >>> sys1.description = "A test for the spin system"
        >>> print(sys1.description)
        A test for the spin system

    transition_pathways: list of :ref:`transition_pathway_api` (optional).
        A list of :ref:`transition_pathway_api` or equivalent dict objects. Each
        transition pathway is a list of :ref:`transition_api` objects. The resulting
        spectrum is a sum of the resonances arising from individual transition pathways.
        The default value is None.

        Example
        -------
        >>> sys1.transition_pathways = [
        ...     [
        ...         {"initial": [-2.5, 0.5], "final": [2.5, 0.5]},
        ...         {"initial": [0.5, 0.5], "final": [-0.5, 0.5]}
        ...     ]
        ... ]
        >>> print(sys1.transition_pathways)
        [|2.5, 0.5⟩⟨-2.5, 0.5| ⟶ |-0.5, 0.5⟩⟨0.5, 0.5|]

        .. note::
            From any given spin system, the list of relevant transition pathways is
            determined by the NMR method. For example, consider a single site
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
                <./../examples/Fitting/plot_2_mrsimFitExample_O17.html>`_
    """

    name: str = None
    label: str = None
    description: str = None
    sites: Union[List[Site], np.ndarray] = []
    couplings: Union[List[Coupling], np.ndarray] = None
    abundance: float = Field(default=100.0, ge=0.0, le=100.0)
    transition_pathways: List = None

    property_unit_types: ClassVar = {"abundance": "dimensionless"}
    property_default_units: ClassVar = {"abundance": "pct"}
    property_units: Dict = {"abundance": "pct"}

    class Config:
        validate_assignment = True
        arbitrary_types_allowed = True

    @validator("transition_pathways")
    def transition_pathways_must_include_transition(cls, v, values):
        if v is None:
            return v
        return [TransitionPathway(item) for item in v]

    @validator("sites")
    def check_sites(cls, v, values):
        if isinstance(v, np.ndarray):
            if not np.all([isinstance(item, Site) for item in v]):
                raise ValueError("All entries must be of type `Site`.")
        return list(v)

    @validator("couplings")
    def check_couplings(cls, v, values):
        if isinstance(v, np.ndarray):
            if not np.all([isinstance(item, Coupling) for item in v]):
                raise ValueError("All entries must be of type `Coupling`.")
        return list(v)

    def get_isotopes(self, spin_I: float = None, symbol: bool = False) -> list:
        """
        An ordered list of :ref:`isotope_api` objects from the sites within the spin
        system corresponding to the given value of spin quantum number `I`. If `I` is
        None, a list of all Isotope objects is returned instead.

        Args:
            float spin_I: An optional spin quantum number. The valid inputs are the
                multiples of 0.5.
            bool symbol: If true, return a list of str with isotope symbols.

        Returns:
            A list of :ref:`isotope_api` objects.

        Example
        -------

        >>> spin_systems.get_isotopes() # three spin systems
        [Isotope(symbol='13C'), Isotope(symbol='1H'), Isotope(symbol='27Al')]
        >>> spin_systems.get_isotopes(symbol=True) # three spin systems
        ['13C', '1H', '27Al']

        >>> spin_systems.get_isotopes(spin_I=0.5) # isotopes with I=0.5
        [Isotope(symbol='13C'), Isotope(symbol='1H')]
        >>> spin_systems.get_isotopes(spin_I=0.5, symbol=True) # isotopes with I=0.5
        ['13C', '1H']

        >>> spin_systems.get_isotopes(spin_I=1.5) # isotopes with I=1.5
        []

        >>> spin_systems.get_isotopes(spin_I=2.5) # isotopes with I=2.5
        [Isotope(symbol='27Al')]
        >>> spin_systems.get_isotopes(spin_I=2.5, symbol=True) # isotopes with I=2.5
        ['27Al']
        """
        isotope_list = allowed_isotopes(spin_I)
        if not symbol:
            return [
                site.isotope
                for site in self.sites
                if site.isotope.symbol in isotope_list
            ]

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
        return get_zeeman_states(self)

    def zeeman_energy_states(self) -> list:
        r"""
        Return a list of all :ref:`zeeman_api` objects of the spin system, where the
        energy states are represented by a list of quantum numbers,

        .. math::
            |\Psi⟩ = [m_1, m_2,.. m_n],

        where :math:`m_i` is the quantum number associated with the :math:`i^\text{th}`
        site within the spin system, and :math:`\Psi` is the energy state.

        Example
        -------

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
        """Returns a list of all possible spin :ref:`transition_api` objects in the
        given spin system.

        Example
        -------

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

        Returns
            A list of :ref:`transition_api` objects.
        """
        transitions = self._all_transitions()
        return TransitionList(
            [
                Transition(initial=item[0].tolist(), final=item[1].tolist())
                for item in transitions
            ]
        )

    @classmethod
    def parse_dict_with_units(cls, py_dict: dict):
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

        if "couplings" in py_dict_copy:
            py_dict_copy["couplings"] = [
                Coupling.parse_dict_with_units(s) for s in py_dict_copy["couplings"]
            ]

        return super().parse_dict_with_units(py_dict_copy)

    # Deprecated
    # def to_freq_dict(self, B0: float) -> dict:
    #     """
    #     Serialize the SpinSystem object to a JSON compliant python dictionary object,
    #     where the attribute value is a numbers expressed in the attribute's default
    #     unit. The default unit for the attributes with respective dimensionalities
    #     are:

    #     - frequency: `Hz`
    #     - angle: `rad`

    #     Args:
    #         float B0: A required macroscopic magnetic flux density in units of T.

    #     Return:
    #         A python dict

    #     Example
    #     -------

    #     >>> pprint(spin_system_1.to_freq_dict(B0=9.4))
    #     {"abundance": 100,
    #      "description": None,
    #      "label": None,
    #      "name": None,
    #      "sites": [{"description": None,
    #                 "isotope": "13C",
    #                 "isotropic_chemical_shift": -2013.1791999999998,
    #                 "label": None,
    #                 "name": None,
    #                 "quadrupolar": None,
    #                 "shielding_antisymmetric": None,
    #                 "shielding_symmetric": {"alpha": None,
    #                                         "beta": None,
    #                                         "eta": 0.5,
    #                                         "gamma": None,
    #                                         "zeta": -1006.5895999999999}}],
    #      "transition_pathways": None}
    #     """
    #     temp_dict = self.dict()
    #     temp_dict["sites"] = [site.to_freq_dict(B0) for site in self.sites]
    #     temp_dict.pop("property_units")
    #     return temp_dict


def allowed_isotopes(spin_I: float = None) -> list:
    """
    List of NMR active isotopes currently supported in ``mrsimulator``.

    Args:
        float spin_I: Optional spin quantum number. The valid values are multiples
            of 0.5. The default is None.

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
