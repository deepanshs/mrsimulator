# -*- coding: utf-8 -*-
import warnings
from copy import deepcopy
from functools import reduce
from itertools import permutations
from typing import ClassVar
from typing import Dict
from typing import List
from typing import Optional
from typing import Union

import csdmpy as cp
import numpy as np
from pydantic import BaseModel
from pydantic import Field
from pydantic import validator

from .abstract_list import TransitionList
from .isotope import Isotope
from .parseable import Parseable
from .post_simulation import PostSimulator
from .transition import Transition


class TransitionQuery(BaseModel):
    """Base TransitionQuery class.

    Attributes:
        P: A list of p symmetry transition, where p = Δm and Δm is the difference
                between spin quantum numbers of the final and initial states.
    """

    P: Optional[Dict] = {"channel-1": [[-1.0]]}
    D: Optional[Dict] = Field(default=None)
    f: Optional[Dict] = Field(default=None)
    transitions: List[Transition] = None

    class Config:
        validate_assignment = True
        arbitrary_types_allowed = True

    def to_dict_with_units(self):
        return {k: v for k, v in self.dict().items() if v is not None}


class Event(Parseable):
    r"""Base Event class defines the spin environment and the transition query for a
    segment of the transition pathway.

    Attributes:
        fraction: A `required` float containing the weight of the frequency
            contribution from the event.
        magnetic_flux_density: An `optional` float containing the macroscopic magnetic
            flux density, :math:`H_0`, of the applied external magnetic field
            during the event in units of T. The default value is ``9.4``.
        rotor_frequency: An `optional` float containing the sample spinning frequency
            :math:`\nu_r`, during the event in units of Hz.
            The default value is ``0``.
        rotor_angle: An `optional` float containing the angle between the
            sample rotation axis and the applied external magnetic field,
            :math:`\theta`, during the event in units of rad.
            The default value is ``0.9553166``, i.e. the magic angle.
        transition_query: An `optional` TransitionQuery object or an equivalent dict
            object listing the queries used in selecting the active transitions
            during the event. Only the active transitions from this query
            contribute to the frequency.
    """

    fraction: float = 1.0
    magnetic_flux_density: Optional[float] = Field(default=9.4, ge=0)
    rotor_frequency: Optional[float] = Field(default=0.0, ge=0)
    # 54.735 degrees = 0.9553166 radians
    rotor_angle: Optional[float] = Field(default=0.9553166, ge=0, le=1.5707963268)
    transition_query: Optional[TransitionQuery] = TransitionQuery()
    user_variables: Optional[List] = None

    property_unit_types: ClassVar = {
        "magnetic_flux_density": "magnetic flux density",
        "rotor_frequency": "frequency",
        "rotor_angle": "angle",
    }

    property_default_units: ClassVar = {
        "magnetic_flux_density": "T",
        "rotor_frequency": "Hz",
        "rotor_angle": "rad",
    }

    property_units: Dict = {
        "magnetic_flux_density": "T",
        "rotor_frequency": "Hz",
        "rotor_angle": "rad",
    }

    class Config:
        validate_assignment = True

    @classmethod
    def parse_dict_with_units(cls, py_dict: dict):
        """
        Parse the physical quantities of an Event object from a python dictionary
        object.

        Args:
            dict py_dict: Dict object
        """
        py_dict_copy = deepcopy(py_dict)
        return super().parse_dict_with_units(py_dict_copy)


class SpectralDimension(Parseable):
    r"""Base SpectralDimension class defines the dimensions of the method.

    Attributes:
        count: An optional integer with the number of points, :math:`N`, along the
            dimension. The default value is 1024.
        spectral_width: An optional float with the spectral width, :math:`\Delta x`,
            along the dimension in units of Hz. The default value is 25 kHz.
        reference_offset: An `optional` float with the reference offset, :math:`x_0`
                along the dimension in units of Hz. The default value is 0.
        origin_offset: An `optional` float with the origin offset (Larmor frequency)
                along the dimension in units of Hz. The default value is None.
        label: An `optional` string label. The default is an empty string.
        events: A `required` list of Event object or an equivalent list of dict objects
                describing the series of events along the spectroscopic dimension.
    """

    count: int = Field(1024, gt=0)
    spectral_width: float = Field(default=25000.0, gt=0)
    reference_offset: Optional[float] = Field(default=0.0)
    origin_offset: Optional[float] = None
    label: Optional[str] = ""
    events: List[Event]

    property_unit_types: ClassVar = {
        "spectral_width": ["frequency", "dimensionless"],
        "reference_offset": ["frequency", "dimensionless"],
        "origin_offset": ["frequency", "dimensionless"],
    }

    property_default_units: ClassVar = {
        "spectral_width": ["Hz", "ppm"],
        "reference_offset": ["Hz", "ppm"],
        "origin_offset": ["Hz", "ppm"],
    }

    property_units: Dict = {
        "spectral_width": "Hz",
        "reference_offset": "Hz",
        "origin_offset": "Hz",
    }

    class Config:
        validate_assignment = True

    @property
    def coordinates_Hz(self) -> np.ndarray:
        r"""
        The grid coordinates along the dimension in units of Hz, evaluated as

        .. math::
            x_\text{Hz} = \left([0, 1, ... N-1] - T\right) \frac{\Delta x}{N} + x_0

        where :math:`T=N/2` and :math:`T=(N-1)/2` for even and odd values of
        :math:`N`, respectively.
        """
        n = self.count
        Tk = int(n / 2)
        increment = self.spectral_width / self.count
        return (np.arange(n) - Tk) * increment + self.reference_offset

    @property
    def coordinates_ppm(self) -> np.ndarray:
        r"""
        The grid coordinates along the dimension as dimension frequency ratio
        in units of ppm. The coordinates are evaluated as

        .. math::
            x_\text{ppm} = \frac{x_\text{Hz}} {x_0 + \omega_0}

        where :math:`\omega_0` is the Larmor frequency.
        """
        if self.origin_offset is None:
            warnings.warn(
                UserWarning(
                    "The coordinates along the dimension without an origin offset "
                    "cannot be converted to dimensionless frequency ratio."
                )
            )
        else:
            denominator = (self.reference_offset + self.origin_offset) / 1e6
            return self.coordinates_Hz / abs(denominator)

    @classmethod
    def parse_dict_with_units(cls, py_dict: dict):
        """
        Parse the physical quantities of a SpectralDimension object from a
        python dictionary object.

        Args:
            dict py_dict: Dict object
        """
        py_dict_copy = deepcopy(py_dict)
        if "events" in py_dict_copy:
            py_dict_copy["events"] = [
                Event.parse_dict_with_units(e) for e in py_dict_copy["events"]
            ]

        return super().parse_dict_with_units(py_dict_copy)

    def to_csdm_dimension(self) -> cp.Dimension:
        """Return the spectral dimension as a CSDM dimension object."""
        increment = self.spectral_width / self.count
        dim = cp.Dimension(
            type="linear",
            count=self.count,
            increment=f"{increment} Hz",
            coordinates_offset=f"{self.reference_offset} Hz",
            complex_fft=True,
            reciprocal={"coordinates_offset": f"{-1/(2*increment)} s"},
        )
        if self.origin_offset is not None:
            dim.origin_offset = f"{self.origin_offset} Hz"
        return dim


class Method(Parseable):
    r"""Base Method class.

    Attributes:
        name: (string) An optional name of the method.
        description: (string) A optional description of the method.
        channels: An required list of isotope symbols over which the given method
            applies, for example, ['1H'].
        spectral_dimensions: An required list of SpectralDimension objects or list of
            equivalent python dictionary objects.
        simulation: A csdm object holding the result of the simulation.
        experiment: A csdm object that optionally holds the experimental data, if
            available.
        post_simulation: An optional dict with post-simulation parameters.
    """
    name: Optional[str] = ""
    description: Optional[str] = ""
    channels: List[str]
    spectral_dimensions: List[SpectralDimension]
    simulation: Optional[Union[cp.CSDM, np.ndarray]]
    post_simulation: Optional[PostSimulator]
    experiment: Optional[Union[cp.CSDM, np.ndarray]]
    # post_simulation: Optional[Dict]

    property_default_units: ClassVar = {
        "magnetic_flux_density": "T",
        "rotor_angle": "rad",
        "rotor_frequency": "Hz",
    }

    property_units: Dict = {
        "magnetic_flux_density": "T",
        "rotor_angle": "rad",
        "rotor_frequency": "Hz",
    }

    class Config:
        validate_assignment = True
        arbitrary_types_allowed = True

    @validator("channels", always=True)
    def validate_channels(cls, v, *, values, **kwargs):
        return [Isotope(symbol=_) for _ in v]

    @validator("experiment", pre=True, always=True)
    def validate_experiment(cls, v, *, values, **kwargs):
        if v is None:
            return None
        if isinstance(v, dict):
            return cp.parse_dict(v)
        if isinstance(v, cp.CSDM):
            return v
        raise ValueError("Unable to read the data.")

    @classmethod
    def parse_dict_with_units(cls, py_dict):
        """
        Parse the physical quantities of the Method object from a python dictionary
        object.

        Args:
            dict py_dict: Dict object
        """
        py_dict_copy = deepcopy(py_dict)
        if "spectral_dimensions" in py_dict_copy:
            py_dict_copy["spectral_dimensions"] = [
                SpectralDimension.parse_dict_with_units(s)
                for s in py_dict_copy["spectral_dimensions"]
            ]
        if "simulation" in py_dict_copy:
            py_dict_copy["simulation"] = cp.parse_dict(py_dict_copy["simulation"])
        if "experiment" in py_dict_copy:
            py_dict_copy["experiment"] = cp.parse_dict(py_dict_copy["experiment"])
        return super().parse_dict_with_units(py_dict_copy)

    def update_spectral_dimension_attributes_from_experiment(self):
        """Update the spectral dimension attributes of the method to match the
        attributes of the experiment from the :attr:`~mrsimulator.Method.experiment`
        attribute."""
        spectral_dims = self.spectral_dimensions
        for i, dim in enumerate(self.experiment.dimensions):
            spectral_dims[i].count = dim.count
            spectral_dims[i].spectral_width = dim.count * dim.increment.to("Hz").value
            spectral_dims[i].reference_offset = dim.coordinates_offset.to("Hz").value
            spectral_dims[i].origin_offset = dim.origin_offset.to("Hz").value

    def to_dict_with_units(self):
        """Parse the class object to a JSON compliant python dictionary object where
        the attribute value with physical quantity is expressed as a string with a
        value and a unit."""
        temp_dict = self.dict(
            exclude={
                "spectral_dimensions",
                "channels",
                "simulation",
                "experiment",
                "property_units",
                "post_simulation",
            }
        )
        temp_dict["spectral_dimensions"] = [
            item.to_dict_with_units() for item in self.spectral_dimensions
        ]
        temp_dict["channels"] = [item.to_dict_with_units() for item in self.channels]
        # if self.simulation is not None:
        #     temp_dict["simulation"] = self.simulation.to_dict(update_timestamp=True)
        # if self.experiment is not None:
        #     temp_dict["experiment"] = self.experiment.to_dict()
        return temp_dict

    def dict(self, **kwargs):
        temp_dict = super().dict(**kwargs)
        if self.post_simulation is not None:
            temp_dict["post_simulation"] = self.post_simulation.dict()
        if self.simulation is not None:
            temp_dict["simulation"] = self.simulation.to_dict(update_timestamp=True)
        if self.experiment is not None:
            temp_dict["experiment"] = self.experiment.to_dict()
        return temp_dict

    def _get_transition_pathways(self, isotopomer):
        selected_transitions = isotopomer._all_transitions()

        segments = []
        for seq in self.spectral_dimensions:
            for ent in seq.events:
                list_of_P = query_permutations(
                    ent.transition_query.to_dict_with_units(),
                    isotope=isotopomer.get_isotopes(),
                    channel=[item.symbol for item in self.channels],
                )
                indexes = P_symmetry_indexes(selected_transitions, list_of_P)
                selected_transitions = selected_transitions[indexes]

                if ent.transition_query.D is not None:
                    list_of_D = query_permutations(
                        ent.transition_query.to_dict_with_units(),
                        isotope=isotopomer.get_isotopes(),
                        channel=[item.symbol for item in self.channels],
                        transition_symmetry="D",
                    )
                    indexes = D_symmetry_indexes(selected_transitions, list_of_D)
                    selected_transitions = selected_transitions[indexes]

                segments += [selected_transitions]
        return segments

    def get_transition_pathways(self, isotopomer) -> np.ndarray:
        """
        Return a list of transition pathways from the given isotopomer that satisfy
        the query selection criterion of the method.

        Args:
            SpinSystem isotopomer: An SpinSystem object.

        Returns: An array of TransitionList objects. Each TransitionList object is a
                transition pathways containing a series of Transition objects.
        """
        segments = self._get_transition_pathways(isotopomer)
        segments = [
            np.asarray(
                TransitionList(
                    [
                        Transition(initial=_[0].tolist(), final=_[1].tolist())
                        for _ in item
                    ]
                )
            )
            for item in segments
        ]
        return cartesian_product(*segments)

    # def get_transition_pathways_old(self, isotopomer):
    #     """
    #     Return a list of transition pathways from the given isotopomer that satisfy
    #     the query criterion of the method.

    #     Args:
    #         isotopomer: An SpinSystem object.
    #     """
    #     transitions = isotopomer.all_transitions
    #     segments = []
    #     for seq in self.spectral_dimensions:
    #         for ent in seq.events:
    #             list_of_P = query_permutations(
    #                 ent.transition_query.to_dict_with_units(),
    #                 isotope=isotopomer.get_isotopes(),
    #                 channel=[item.symbol for item in self.channels],
    #             )
    #             P_segment = []
    #             # delta_P = transitions[:,:]
    #             # for symmetry in list_of_P:
    #             #     P_segment += transitions.filter(P=symmetry)

    #             # if ent.transition_query.D != None:
    #             #     list_of_D = query_permutations(
    #             #         ent.transition_query.to_dict_with_units(),
    #             #         isotope=isotopomer.get_isotopes(),
    #             #         channel=[item.symbol for item in self.channels],
    #             #         transition_symmetry = "D"
    #             #     )

    #             #     D_segment = []
    #             #     for D_symmetry in list_of_D:
    #             #         D_segment += transitions.filter(D = D_symmetry)
    #             #     print('list of D: ', list_of_D)
    #             #     print('D_segment: ', D_segment)

    #             for symmetry in list_of_P:
    #                 if ent.transition_query.D is None:
    #                     P_segment += transitions.filter(P=symmetry)
    #                 elif ent.transition_query.D is not None:
    #                     list_of_D = query_permutations(
    #                         ent.transition_query.to_dict_with_units(),
    #                         isotope=isotopomer.get_isotopes(),
    #                         channel=[item.symbol for item in self.channels],
    #                         transition_symmetry="D",
    #                     )
    #                     D_transition = []
    #                     [
    #                         D_transition.append(x)
    #                         for x in list_of_D
    #                         if x not in D_transition
    #                     ]
    #                     # D_segment = []
    #                     for D_symmetry in D_transition:
    #                         P_segment += transitions.filter(P=symmetry, D=D_symmetry)
    #                     # print('list of D: ', D_transition)
    #                     # for symmetry in list_of_D:
    #                     #     D_segment += transitions.filter(D=symmetry)
    #                     # print('D_segment: ', D_segment)

    #             segments.append(np.asarray(P_segment))  # append the intersection
    #     return cartesian_product(*segments)

    def apodize(self, **kwargs):
        """Returns the result of passing the selected apodization function .

        Args:
            csdm: simulation object

        Returns:
            A Numpy array
        """

        csdm = self.simulation
        for dim in csdm.dimensions:
            dim.to("Hz", "nmr_frequency_ratio")
        apo = self.post_simulation.apodization

        sum_ = 0

        for item in apo:
            sum_ += item.apodize(csdm)

        for dim in csdm.dimensions:
            dim.to("ppm", "nmr_frequency_ratio")

        return self.post_simulation.scale * sum_


def cartesian_product(*arrays):
    la = len(arrays)
    dtype = np.result_type(*arrays)
    arr = np.empty([len(a) for a in arrays] + [la], dtype=dtype)
    for i, a in enumerate(np.ix_(*arrays)):
        arr[..., i] = a
    return arr.reshape(-1, la)


def P_symmetry_indexes(transitions, list_of_P):
    P = transitions[:, 1, :] - transitions[:, 0, :]
    n = P.shape[1]
    return reduce(
        np.union1d,
        [
            reduce(np.intersect1d, [np.where(P[:, i] == search[i]) for i in range(n)])
            for search in list_of_P
        ],
        np.asarray([], dtype=np.int64),
    )


def D_symmetry_indexes(transitions, list_of_D):
    D = transitions[:, 1, :] ** 2 - transitions[:, 0, :] ** 2
    return reduce(
        np.union1d,
        [
            reduce(
                np.intersect1d,
                [np.where(D[:, i] == search[i]) for i in range(D.shape[1])],
            )
            for search in list_of_D
        ],
        np.asarray([], dtype=np.int64),
    )


def get_iso_dict(channel, isotope):
    """
        Parse the isotopomer sites to determine indices of each isotope that
        is part of the method channel.

        Args:
            channel: List object
            isotope: List object

    """
    iso_dict = {}

    # determine channels for P
    for i, item in enumerate(isotope):
        if item in channel and item not in iso_dict:
            iso_dict[item] = [i]
        elif item in iso_dict:
            iso_dict[item].append(i)

    return iso_dict


def query_permutations(query, isotope, channel, transition_symmetry="P"):
    """
        Determines the transition symmetries that are involved in a given
        transition query.

        Args:
            query: Dict object
            channel: List object
            isotope: List object
            transition_symmetry: str object. Derived from a transition query

    """

    P_permutated = []
    iso_dict = get_iso_dict(channel=channel, isotope=isotope)
    query_short = query[transition_symmetry]
    for i, items in enumerate(query_short):
        # Check if method isotope is in the isotopomer
        if channel[i] not in iso_dict:
            print(
                f"Method/channel isotope mismatch. Channel asks for {channel[i]} "
                f"but is not in {isotope}"
            )
            return []

        temp_P = []
        for k in range(len(query_short[items])):
            # Check transition query doesn't require more isotopes than present
            if len(query_short[items][k]) > len(iso_dict[channel[i]]):
                print("Failed: Transition query larger than channel")
                return []
            elif len(query_short[items][k]) <= len(iso_dict[channel[i]]):
                query_short[items][k] += (
                    len(iso_dict[channel[i]]) - len(query_short[items][k])
                ) * [k]
            temp_P += list(permutations(query_short[items][k]))
        P_permutated += [temp_P]
        # P_permutated += [list(permutations(query_short[items][k]))]

    transition_symmetry_from_query = []
    for i, iso_trans_symmetry in enumerate(P_permutated):
        # creating transition symmetries isotope by isotope
        temp_transitions = []
        for transition in iso_trans_symmetry:
            P_expanded = len(isotope) * [0]
            for j, item in enumerate(transition):
                # filling indices of isotopomer with a sites transition symmetries
                P_expanded[iso_dict[channel[i]][j]] = item

            if transition_symmetry_from_query == []:
                temp_transitions.append(P_expanded)
            else:
                # Each isotope is added to the previous isotope to create the
                # full transition symmetry
                for k, intermediate in enumerate(transition_symmetry_from_query[-1]):
                    temp_transitions.append(
                        [sum(x) for x in zip(intermediate, P_expanded)]
                    )

        transition_symmetry_from_query.append(temp_transitions)

    return transition_symmetry_from_query[-1]
