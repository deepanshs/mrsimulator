# -*- coding: utf-8 -*-
from copy import deepcopy
from typing import ClassVar
from typing import Dict
from typing import List
from typing import Union

import csdmpy as cp
import numpy as np
from mrsimulator.spin_system.isotope import Isotope
from mrsimulator.transition import Transition
from mrsimulator.transition.transition_list import TransitionList
from mrsimulator.utils.parseable import Parseable
from pydantic import validator

from .spectral_dimension import SpectralDimension
from .utils import cartesian_product
from .utils import D_symmetry_indexes
from .utils import P_symmetry_indexes
from .utils import query_permutations

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


class Method(Parseable):
    r"""Base Method class. A method class represents the NMR method.

    Attributes
    ----------

    channels: list (optional).
        The value is a list of isotope symbols over which the given method applies.
        An isotope symbol is given as a string with the atomic number followed by its
        atomic symbol, for example, '1H', '13C', and '33S'. The default is an empty
        list.
        The number of isotopes in a `channel` depends on the method. For example, a
        `BlochDecaySpectrum` method is a single channel method, in which case, the
        value of this attribute is a list with a single isotope symbol, ['13C'].

        Example
        -------

        >>> bloch = Method()
        >>> bloch.channels = ['1H']

    spectral_dimensions: list of :ref:`spectral_dim_api` or dict objects (optional).
        The number of spectral dimensions depends on the given method. For example, a
        `BlochDecaySpectrum` method is a one-dimensional method and thus requires a
        single spectral dimension. The default is a single default
        :ref:`spectral_dim_api` object.

        Example
        -------

        >>> bloch = Method()
        >>> bloch.spectral_dimensions = [SpectralDimension(count=8, spectral_width=50)]
        >>> # or equivalently
        >>> bloch.spectral_dimensions = [{'count': 8, 'spectral_width': 50}]

    simulation: CSDM or ndarray (N/A).
        An object holding the result of the simulation. The initial value of this
        attribute is None. A value is assigned to this attribute when you run the
        simulation using the :meth:`~mrsimulator.Simulator.run` method.

    experiment: CSDM or ndarray (optional).
        An object holding the experimental measurement for the given method, if
        available. The default value is None.

        Example
        -------

        >>> bloch.experiment = my_data # doctest: +SKIP

    name: str (optional).
        The value is the name or id of the method. The default value is None.

        Example
        -------

        >>> bloch.name = 'BlochDecaySpectrum'
        >>> bloch.name
        'BlochDecaySpectrum'

    label: str (optional).
        The value is a label for the method. The default value is None.

        Example
        -------

        >>> bloch.label = 'One pulse acquired spectrum'
        >>> bloch.label
        'One pulse acquired spectrum'

    description: str (optional).
        The value is a description of the method. The default value is None.

        Example
        -------

        >>> bloch.description = 'Huh!'
        >>> bloch.description
        'Huh!'

    """
    name: str = None
    label: str = None
    description: str = None
    channels: List[str] = []
    spectral_dimensions: List[SpectralDimension] = [{}]
    simulation: Union[cp.CSDM, np.ndarray] = None
    experiment: Union[cp.CSDM, np.ndarray] = None

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
        Parse the physical quantity from a dictionary representation of the Method
        object, where the physical quantity is expressed as a string with a number and
        a unit.

        Args:
            dict py_dict: A python dict representation of the Method object.

        Returns:
            A :ref:`method_api` object.
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
        """
        Parse the class object to a JSON compliant python dictionary object where
        the attribute value with physical quantity is expressed as a string with a
        value and a unit.

        Returns:
            A python dict object.
        """
        temp_dict = {}
        items = ["name", "label", "description"]
        for en in items:
            value = self.__getattribute__(en)
            if value is not None:
                temp_dict[en] = self.__getattribute__(en)

        temp_dict["spectral_dimensions"] = [
            item.to_dict_with_units() for item in self.spectral_dimensions
        ]
        temp_dict["channels"] = [item.to_dict_with_units() for item in self.channels]
        if self.simulation is not None:
            temp_dict["simulation"] = self.simulation.to_dict(update_timestamp=True)
        if self.experiment is not None:
            temp_dict["experiment"] = self.experiment.to_dict()
        return temp_dict

    def dict(self, **kwargs):
        temp_dict = super().dict(**kwargs)
        if self.simulation is not None:
            temp_dict["simulation"] = self.simulation.to_dict(update_timestamp=True)
        if self.experiment is not None and isinstance(self.experiment, cp.CSDM):
            temp_dict["experiment"] = self.experiment.to_dict()
        return temp_dict

    def _get_transition_pathways(self, spin_system):
        all_transitions = spin_system._all_transitions()

        segments = []
        for seq in self.spectral_dimensions:
            selected_transitions = all_transitions[:]
            for ent in seq.events:
                list_of_P = query_permutations(
                    ent.transition_query.to_dict_with_units(),
                    isotope=spin_system.get_isotopes(),
                    channel=[item.symbol for item in self.channels],
                )
                indexes = P_symmetry_indexes(selected_transitions, list_of_P)
                selected_transitions = selected_transitions[indexes]

                if ent.transition_query.D is not None:
                    list_of_D = query_permutations(
                        ent.transition_query.to_dict_with_units(),
                        isotope=spin_system.get_isotopes(),
                        channel=[item.symbol for item in self.channels],
                        transition_symmetry="D",
                    )
                    indexes = D_symmetry_indexes(selected_transitions, list_of_D)
                    selected_transitions = selected_transitions[indexes]

                segments += [selected_transitions]
        return segments

    def _get_transition_pathways_np(self, spin_system):
        segments = self._get_transition_pathways(spin_system)
        segments_index = [np.arange(item.shape[0]) for item in segments]
        cartesian_index = cartesian_product(*segments_index)
        return [
            [segments[i][j] for i, j in enumerate(item)] for item in cartesian_index
        ]

    def get_transition_pathways(self, spin_system) -> np.ndarray:
        """
        Return a list of transition pathways from the given spin system that satisfy
        the query selection criterion of the method.

        Args:
            SpinSystem spin_system: A SpinSystem object.

        Returns:
            An array of TransitionList objects. Each TransitionList object is a
            transition pathways containing a series of Transition objects.
        """
        segments = self._get_transition_pathways_np(spin_system)
        return np.asarray(
            [
                TransitionList(
                    [
                        Transition(initial=tr[0].tolist(), final=tr[1].tolist())
                        for tr in item
                    ]
                )
                for item in segments
            ]
        )
