# -*- coding: utf-8 -*-
from copy import deepcopy
from os import path
from typing import ClassVar
from typing import Dict
from typing import List
from typing import Union

import csdmpy as cp
import numpy as np
from monty.serialization import loadfn
from mrsimulator.spin_system.isotope import Isotope
from mrsimulator.transition import Transition
from mrsimulator.transition import TransitionPathway
from mrsimulator.utils.parseable import Parseable
from pydantic import validator

from .event import Event
from .spectral_dimension import SpectralDimension
from .utils import cartesian_product
from .utils import D_symmetry_indexes
from .utils import expand_spectral_dimension_object
from .utils import P_symmetry_indexes
from .utils import query_permutations

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"

MODULE_DIR = path.dirname(path.abspath(__file__))
NAMED_METHODS = loadfn(path.join(MODULE_DIR, "named_methods.json"))["named_methods"]


class Method(Parseable):
    r"""Base Method class. A method class represents the NMR method.

    Attributes
    ----------

    channels: List[str] (optional).
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

    spectral_dimensions: List[:ref:`spectral_dim_api`] or List[dict] (optional).
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
        Name or id of the method. The default value is None.

        Example
        -------

        >>> bloch.name = 'BlochDecaySpectrum'
        >>> bloch.name
        'BlochDecaySpectrum'

    label: str (optional).
        Label for the method. The default value is None.

        Example
        -------

        >>> bloch.label = 'One pulse acquired spectrum'
        >>> bloch.label
        'One pulse acquired spectrum'

    description: str (optional).
        A description of the method. The default value is None.

        Example
        -------

        >>> bloch.description = 'Huh!'
        >>> bloch.description
        'Huh!'

    affine_matrix: np.ndarray or 2D list (optional)
        A (`n` x `n`) affine transformation matrix, where `n` is the number of
        spectral_dimensions. If provided, the corresponding affine transformation is
        applied to the computed frequencies. The default is None, i.e., no
        transformation is applied.

        Example
        -------

        >>> method = Method2D()
        >>> method.affine_matrix = [[1, -1], [0, 1]]
        >>> print(method.affine_matrix)
        [[ 1 -1]
         [ 0  1]]
    """
    name: str = None
    label: str = None
    description: str = None
    channels: List[str] = []
    spectral_dimensions: List[SpectralDimension] = [SpectralDimension()]
    affine_matrix: Union[np.ndarray, List] = None
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

    def __eq__(self, other):
        if not isinstance(other, Method):
            return False
        check = [
            self.name == other.name,
            self.label == other.label,
            self.description == other.description,
            self.channels == other.channels,
            self.spectral_dimensions == other.spectral_dimensions,
            np.all(self.affine_matrix == other.affine_matrix),
            self.simulation == other.simulation,
            self.experiment == other.experiment,
        ]
        if np.all(check):
            return True
        return False

    @validator("channels", always=True)
    def validate_channels(cls, v, *, values, **kwargs):
        return [Isotope(symbol=_) for _ in v]

    @staticmethod
    def __check_csdm__(data):
        if data is None:
            return None
        if isinstance(data, dict):
            return cp.parse_dict(data)
        if isinstance(data, cp.CSDM):
            return data
        raise ValueError("Unable to read the data.")

    @validator("experiment", pre=True, always=True)
    def validate_experiment(cls, v, *, values, **kwargs):
        return cls.__check_csdm__(v)

    @validator("simulation", pre=True, always=True)
    def validate_simulation(cls, v, *, values, **kwargs):
        if isinstance(v, np.ndarray):
            return v
        return cls.__check_csdm__(v)

    @validator("affine_matrix", pre=True, always=True)
    def validate_affine_matrix(cls, v, *, values, **kwargs):
        if v is None:
            return None
        v = np.asarray(v)
        dim_len = len(values["spectral_dimensions"])
        if v.size != dim_len ** 2:
            raise ValueError(f"Expecting a {dim_len}x{dim_len} affine matrix.")
        if v.ravel()[0] == 0:
            raise ValueError("The first element of the affine matrix cannot be zero.")
        return v

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
            py_dict_copy = expand_spectral_dimension_object(py_dict_copy)
            py_dict_copy["spectral_dimensions"] = [
                SpectralDimension.parse_dict_with_units(s)
                for s in py_dict_copy["spectral_dimensions"]
            ]

        if "simulation" in py_dict_copy:
            if py_dict_copy["simulation"] is not None:
                py_dict_copy["simulation"] = cp.parse_dict(py_dict_copy["simulation"])
        if "experiment" in py_dict_copy:
            if py_dict_copy["experiment"] is not None:
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

    def dict(self, **kwargs):
        mth = super().dict(**kwargs)
        if self.simulation is not None:
            mth["simulation"] = self.simulation.to_dict(update_timestamp=True)
        if self.experiment is not None and isinstance(self.experiment, cp.CSDM):
            mth["experiment"] = self.experiment.to_dict()
        return mth

    def json(self) -> dict:
        """Parse the class object to a JSON compliant python dictionary object, where
        the attribute value with physical quantity is expressed as a string with a
        value and a unit.

        Returns:
            A python dict object.
        """

        mth = {_: self.__getattribute__(_) for _ in ["name", "label", "description"]}
        mth["channels"] = [item.json() for item in self.channels]
        mth["spectral_dimensions"] = [item.json() for item in self.spectral_dimensions]

        # add global parameters
        ev0 = self.spectral_dimensions[0].events[0]
        evt_d = Event.property_default_units
        global_ = {k: f"{ev0.__getattribute__(k)} {u}" for k, u in evt_d.items()}
        mth.update(global_)

        named = True if mth["name"] in NAMED_METHODS else False
        for dim in mth["spectral_dimensions"]:
            for ev in dim["events"]:
                # remove event objects with global values.
                _ = [ev.pop(k) if ev[k] == v else 0 for k, v in global_.items()]

                # remove transition query objects for named methods
                _ = ev.pop("transition_query") if named else 0

            if dim["events"] == [{} for _ in dim["events"]]:
                dim.pop("events")

        afm = self.affine_matrix
        mth["affine_matrix"] = None if afm is None else afm.tolist()

        sim = self.simulation
        mth["simulation"] = None if sim is None else sim.to_dict(update_timestamp=True)

        exp = self.experiment
        mth["experiment"] = None if exp is None else exp.to_dict()

        _ = [mth.pop(item) for item in [k for k, v in mth.items() if v is None]]
        return mth

    # def _simplify_events_json(self, py_dict):
    #     named = True if py_dict["name"] in NAMED_METHODS else False
    #     for dim in py_dict["spectral_dimensions"]:
    #         for ev in dim["events"]:
    #             # remove event objects with global values.
    #             for key, val in zip(list_g, global_):
    #                 _ = ev.pop(key) if ev[key] == val else 0

    #             # remove transition query objects for named methods
    #             _ = ev.pop("transition_query") if named else 0

    #         if dim["events"] == [{} for _ in range(len(dim["events"]))]:
    #             dim.pop("events")

    # def _get_symmetry_pathways(self, spin_system):
    #     list_of_P = []
    #     list_of_D = []
    #     for dim in self.spectral_dimensions:
    #         for ent in dim.events:
    #             list_of_P.append(
    #                 query_permutations(
    #                     ent.transition_query.dict(),
    #                     isotope=spin_system.get_isotopes(symbol=True),
    #                     channel=[item.symbol for item in self.channels],
    #                 )
    #             )
    #             if ent.transition_query.D is not None:
    #                 list_of_D.append(
    #                     query_permutations(
    #                         ent.transition_query.dict(),
    #                         isotope=spin_system.get_isotopes(symbol=True),
    #                         channel=[item.symbol for item in self.channels],
    #                         transition_symmetry="D",
    #                     )
    #                 )

    #     return {"P": list_of_P, "D": list_of_D}

    def _get_transition_pathways(self, spin_system):
        all_transitions = spin_system._all_transitions()

        segments = []
        for dim in self.spectral_dimensions:
            for ent in dim.events:
                # query the transitions for P symmetry
                selected_transitions = all_transitions[:]
                list_of_P = query_permutations(
                    ent.transition_query.dict(),
                    isotope=spin_system.get_isotopes(symbol=True),
                    channel=[item.symbol for item in self.channels],
                )
                indexes = P_symmetry_indexes(selected_transitions, list_of_P)
                selected_transitions = selected_transitions[indexes]

                # query the transitions for D symmetry
                if ent.transition_query.D is not None:
                    list_of_D = query_permutations(
                        ent.transition_query.dict(),
                        isotope=spin_system.get_isotopes(symbol=True),
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

    def get_transition_pathways(self, spin_system) -> List[TransitionPathway]:
        """
        Return a list of transition pathways from the given spin system that satisfy
        the query selection criterion of the method.

        Args:
            SpinSystem spin_system: A SpinSystem object.

        Returns:
            An array of :ref:`transition_pathway_api` objects. Each TransitionPathway
            object is an ordered collection of Transition objects.

        Example:
            >>> from mrsimulator import SpinSystem
            >>> from mrsimulator.methods import ThreeQ_VAS
            >>> sys = SpinSystem(sites=[{'isotope': '27Al'}, {'isotope': '29Si'}])
            >>> method = ThreeQ_VAS(channels=['27Al'])
            >>> pprint(method.get_transition_pathways(sys))
            [|1.5, -0.5⟩⟨-1.5, -0.5| ⟶ |-0.5, -0.5⟩⟨0.5, -0.5|,
             |1.5, -0.5⟩⟨-1.5, -0.5| ⟶ |-0.5, 0.5⟩⟨0.5, 0.5|,
             |1.5, 0.5⟩⟨-1.5, 0.5| ⟶ |-0.5, -0.5⟩⟨0.5, -0.5|,
             |1.5, 0.5⟩⟨-1.5, 0.5| ⟶ |-0.5, 0.5⟩⟨0.5, 0.5|]
        """
        segments = self._get_transition_pathways_np(spin_system)
        return [
            TransitionPathway(
                [
                    Transition(initial=tr[0].tolist(), final=tr[1].tolist())
                    for tr in item
                ]
            )
            for item in segments
        ]

    def shape(self) -> tuple:
        """The shape of the method's spectral dimension array.

        Returns:
            tuple

        Example:
            >>> from mrsimulator.methods import Method2D
            >>> method = Method2D(spectral_dimensions=[{'count': 40}, {'count': 10}])
            >>> method.shape()
            (40, 10)
        """
        return tuple([item.count for item in self.spectral_dimensions])
