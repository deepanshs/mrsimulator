from copy import deepcopy
from typing import ClassVar
from typing import Dict
from typing import List
from typing import Union

import csdmpy as cp
import matplotlib as mpl
import numpy as np
import pandas as pd
from mrsimulator.base_model import transition_connect_factor
from mrsimulator.spin_system.isotope import Isotope
from mrsimulator.transition import SymmetryPathway
from mrsimulator.transition import Transition
from mrsimulator.transition import TransitionPathway
from mrsimulator.utils.parseable import Parseable
from pydantic import Field
from pydantic import PrivateAttr
from pydantic import validator

from .event import MixingEvent  # noqa: F401
from .event import SpectralEvent  # noqa: F401
from .plot import plot as _plot
from .spectral_dimension import CHANNELS
from .spectral_dimension import SpectralDimension
from .utils import cartesian_product
from .utils import check_for_at_least_one_event
from .utils import check_for_number_of_spectral_dimensions
from .utils import check_spectral_dimensions_are_dict
from .utils import mixing_query_connect_map
from .utils import to_euler_list

# from .utils import convert_transition_query

# from .event import ConstantDurationEvent  # noqa: F401


__author__ = ["Deepansh J. Srivastava", "Matthew D. Giammar"]
__email__ = ["srivastava.89@osu.edu", "giammar.7@buckeyemail.osu.edu"]


class Method(Parseable):
    r"""Base Method class. A method class represents the NMR method.

    Attributes
    ----------

    channels:
        The value is a list of isotope symbols over which the given method applies.
        An isotope symbol is given as a string with the atomic number followed by its
        atomic symbol, for example, '1H', '13C', and '33S'. The default is an empty
        list.
        The number of isotopes in a `channel` depends on the method. For example, a
        `BlochDecaySpectrum` method is a single channel method, in which case, the
        value of this attribute is a list with a single isotope symbol, ['13C'].

        Example
        -------

        >>> bloch = Method(channels=['1H'], spectral_dimensions=[{}])
        >>> bloch.channels = ['13C']  # Change channels

    spectral_dimensions:
        The number of spectral dimensions depends on the given method. For example, a
        `BlochDecaySpectrum` method is a one-dimensional method and thus requires a
        single spectral dimension.

        Example
        -------

        >>> bloch = Method(channels=['1H'], spectral_dimensions=[
        ...     SpectralDimension(count=8, spectral_width=50)
        ... ])
        >>> # or equivalently
        >>> bloch = Method(channels=['1H'], spectral_dimensions=[
        ...     {"count": 8, "spectral_width": 50}
        ... ])

    simulation:
        An object holding the result of the simulation. The initial value of this
        attribute is None. A value is assigned to this attribute when you run the
        simulation using the :py:meth:`~mrsimulator.Simulator.run` method.

    experiment:
        An object holding the experimental measurement for the given method, if
        available. The default value is None.

        Example
        -------

        >>> bloch.experiment = my_dataset # doctest: +SKIP

    name:
        Name or id of the method. The default value is None.

        Example
        -------

        >>> bloch.name = 'BlochDecaySpectrum'
        >>> bloch.name
        'BlochDecaySpectrum'

    label:
        Label for the method. The default value is None.

        Example
        -------

        >>> bloch.label = 'One pulse acquired spectrum'
        >>> bloch.label
        'One pulse acquired spectrum'

    description:
        A description of the method. The default value is None.

        Example
        -------

        >>> bloch.description = 'Huh!'
        >>> bloch.description
        'Huh!'

    affine_matrix:
        A (`n` x `n`) affine transformation matrix, where `n` is the number of
        spectral_dimensions. If provided, the corresponding affine transformation is
        applied to the computed frequencies. The default is None, i.e., no
        transformation is applied.

        Example
        -------

        >>> method = Method(channels=['1H'], spectral_dimensions=[{}, {}]) # 2D method
        >>> method.affine_matrix = [[1, -1], [0, 1]]
        >>> print(method.affine_matrix)
        [[1, -1], [0, 1]]
    """
    channels: List[str]
    spectral_dimensions: List[SpectralDimension] = [SpectralDimension()]
    affine_matrix: List = None
    simulation: Union[cp.CSDM, np.ndarray] = None
    experiment: Union[cp.CSDM, np.ndarray] = None

    # global attributes
    magnetic_flux_density: float = Field(default=9.4, ge=0.0)
    rotor_frequency: float = Field(default=0.0, ge=0.0)
    rotor_angle: float = Field(default=0.9553166181245, ge=0.0, le=1.5707963268)

    # private attributes
    _named_method: bool = PrivateAttr(False)
    _metadata: dict = PrivateAttr({})

    property_unit_types: ClassVar[Dict] = {
        "magnetic_flux_density": "magnetic flux density",
        "rotor_frequency": "frequency",
        "rotor_angle": "angle",
    }

    property_default_units: ClassVar[Dict] = {
        "magnetic_flux_density": "T",
        "rotor_angle": "rad",
        "rotor_frequency": "Hz",
    }

    property_units: Dict = {
        "magnetic_flux_density": "T",
        "rotor_angle": "rad",
        "rotor_frequency": "Hz",
    }
    test_vars: ClassVar[Dict] = {"channels": ["1H"]}

    class Config:
        extra = "forbid"
        validate_assignment = True
        arbitrary_types_allowed = True

    @validator("channels", always=True)
    def validate_channels(cls, v, *, values, **kwargs):
        return [Isotope(symbol=_) for _ in v]

    def __init__(self, **kwargs):
        Method.check(kwargs)
        super().__init__(**kwargs)
        _ = [
            setattr(ev, item, getattr(self, item))
            for sd in self.spectral_dimensions
            for ev in sd.events
            for item in self.property_units.keys()
            if hasattr(ev, item) and getattr(ev, item) is None
        ]
        # Check for only one event with 0 < rotor_freq < infinity
        speeds = [
            ev.rotor_frequency
            for sd in self.spectral_dimensions
            for ev in sd.events
            if ev.__class__.__name__ not in ["MixingEvent", "ConstantTimeEvent"]
        ]
        speeds = [sp for sp in speeds if 0 < sp < 1e12]
        if len(speeds) > 1:
            raise NotImplementedError(
                f"Sideband-sideband correlation is not yet supported in mrsimulator. "
                f"Only one event with non-zero and finite `rotor_frequency` is allowed "
                f"in a method. Found {len(speeds)} events with rotor frequencies "
                f"{speeds}."
            )

    @classmethod
    def check(cls, kwargs, is_named_method=False, ndim=None):
        check_for_number_of_spectral_dimensions(kwargs, is_named_method, ndim)
        sd_is_dict = check_spectral_dimensions_are_dict(kwargs)
        if sd_is_dict and not is_named_method:
            check_for_at_least_one_event(kwargs)

    @staticmethod
    def __check_csdm__(dataset):
        if dataset is None:
            return None
        if isinstance(dataset, dict):
            return cp.parse_dict(dataset)
        if isinstance(dataset, cp.CSDM):
            return dataset
        raise ValueError("Unable to read the dataset.")

    @validator("experiment", pre=True, always=True)
    def validate_experiment(cls, v, *, values, **kwargs):
        return cls.__check_csdm__(v)

    @validator("simulation", pre=True, always=True)
    def validate_simulation(cls, v, *, values, **kwargs):
        if isinstance(v, np.ndarray):
            return v
        return cls.__check_csdm__(v)

    @validator("affine_matrix", always=True)
    def validate_affine_matrix(cls, v, *, values, **kwargs):
        if v is None:
            return
        v1 = np.asarray(v)
        dim_len = len(values["spectral_dimensions"])
        if v1.size != dim_len**2:
            raise ValueError(f"Expecting a {dim_len}x{dim_len} affine matrix.")
        if v1.ravel()[0] == 0:
            raise ValueError("The first element of the affine matrix cannot be zero.")
        return v

    @classmethod
    def parse_dict_with_units(cls, py_dict):
        """Parse the physical quantity from a dictionary representation of the Method
        object, where the physical quantity is expressed as a string with a number and
        a unit.

        Args:
            dict py_dict: A python dict representation of the Method object.

        Returns:
            A :ref:`method_api` object.
        """
        py_dict_copy = deepcopy(py_dict)

        if "spectral_dimensions" in py_dict_copy:
            Method.expand_spectral_dimension_object(py_dict_copy)
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

    @staticmethod
    def expand_spectral_dimension_object(py_dict):
        glb = {}
        _ = [
            glb.update({item: py_dict[item]})
            for item in Method.property_unit_types.keys()
            if item in py_dict.keys()
        ]
        glb_keys = set(glb.keys())

        _ = [
            (
                None if "events" in dim else dim.update({"events": [{}]}),
                [
                    ev.update({k: glb[k]})
                    for ev in dim["events"]
                    for k in glb
                    if k not in set(ev.keys()).intersection(glb_keys)
                ],
            )
            for dim in py_dict["spectral_dimensions"]
        ]

    def dict(self, **kwargs):
        mth = super().dict(**kwargs)
        if isinstance(self.simulation, cp.CSDM):
            mth["simulation"] = self.simulation.to_dict(update_timestamp=True)
        if isinstance(self.experiment, cp.CSDM):
            mth["experiment"] = self.experiment.to_dict()
        return mth

    def json(self, units=True) -> dict:
        """Parse the class object to a JSON compliant python dictionary object.

        Args:
            units: If true, the attribute value is a physical quantity expressed as a
                string with a number and a unit, else a float.

        Returns: dict
        """
        mth = {k: getattr(self, k) for k in ["name", "label", "description"]}
        mth["channels"] = [item.json() for item in self.channels]
        mth["spectral_dimensions"] = [
            item.json(units=units) for item in self.spectral_dimensions
        ]

        # add global parameters
        global_ = (
            {k: f"{getattr(self, k)} {u}" for k, u in self.property_units.items()}
            if units
            else {k: getattr(self, k) for k in self.property_units}
        )
        mth.update(global_)

        # remove event objects with global values.
        _ = [
            [ev.pop(k) if k in ev and ev[k] == v else 0 for k, v in global_.items()]
            for dim in mth["spectral_dimensions"]
            for ev in dim["events"]
        ]

        # if self._named_method:
        #     _ = [dim.pop("events") for dim in mth["spectral_dimensions"]]

        mth["affine_matrix"] = self.affine_matrix

        sim = self.simulation
        mth["simulation"] = None if sim is None else sim.to_dict(update_timestamp=True)

        exp = self.experiment
        mth["experiment"] = None if exp is None else exp.to_dict()

        _ = [mth.pop(item) for item in [k for k, v in mth.items() if v is None]]
        return mth

    def get_symmetry_pathways(self, symmetry_element: str) -> List[SymmetryPathway]:
        """Return a list of symmetry pathways of the method.

        Args:
            str symmetry_element: The  symmetry element, 'P' or 'D'.

        Returns:
            A list of :ref:`symmetry_pathway_api` objects.

        **Single channel example**

        Example:
            >>> from mrsimulator.method import Method
            >>> method = Method(
            ...     channels=['1H'],
            ...     spectral_dimensions=[
            ...         {
            ...             "events": [
            ...                 {
            ...                     "fraction": 0.5,
            ...                     "transition_queries": [{"ch1": {"P": [1]}}]
            ...                 },
            ...                 {
            ...                     "fraction": 0.5,
            ...                     "transition_queries": [{"ch1": {"P": [0]}}]
            ...                 }
            ...             ],
            ...         },
            ...         {
            ...             "events": [
            ...                 {"transition_queries": [{"ch1": {"P": [-1]}}]},
            ...             ],
            ...         }
            ...     ]
            ... )
            >>> pprint(method.get_symmetry_pathways("P"))
            [SymmetryPathway(
                ch1(1H): [1] ⟶ [0] ⟶ [-1]
                total: 1.0 ⟶ 0.0 ⟶ -1.0
            )]

        **Dual channels example**

        Example:
            >>> from mrsimulator.method import Method
            >>> method = Method(
            ...     channels=['1H', '13C'],
            ...     spectral_dimensions=[
            ...         {
            ...             "events": [{
            ...                 "fraction": 0.5,
            ...                 "transition_queries": [
            ...                     {"ch1": {"P": [1]}},
            ...                     {"ch1": {"P": [-1]}},
            ...                 ]
            ...             },
            ...             {
            ...                 "fraction": 0.5,
            ...                 "transition_queries": [  # selecting double quantum
            ...                     {"ch1": {"P": [-1]}, "ch2": {"P": [-1]}},
            ...                     {"ch1": {"P": [1]}, "ch2": {"P": [1]}},
            ...                 ]
            ...             }],
            ...         },
            ...         {
            ...             "events": [{
            ...                 "transition_queries": [ # selecting single quantum
            ...                     {"ch1": {"P": [-1]}},
            ...                 ]
            ...             }],
            ...         }
            ...     ]
            ... )
            >>> pprint(method.get_symmetry_pathways("P"))
            [SymmetryPathway(
                ch1(1H): [1] ⟶ [-1] ⟶ [-1]
                ch2(13C): None ⟶ [-1] ⟶ None
                total: 1.0 ⟶ -2.0 ⟶ -1.0
            ),
             SymmetryPathway(
                ch1(1H): [1] ⟶ [1] ⟶ [-1]
                ch2(13C): None ⟶ [1] ⟶ None
                total: 1.0 ⟶ 2.0 ⟶ -1.0
            ),
             SymmetryPathway(
                ch1(1H): [-1] ⟶ [-1] ⟶ [-1]
                ch2(13C): None ⟶ [-1] ⟶ None
                total: -1.0 ⟶ -2.0 ⟶ -1.0
            ),
             SymmetryPathway(
                ch1(1H): [-1] ⟶ [1] ⟶ [-1]
                ch2(13C): None ⟶ [1] ⟶ None
                total: -1.0 ⟶ 2.0 ⟶ -1.0
            )]
        """
        sym_path = [
            dim._get_symmetry_pathways(symmetry_element)
            for dim in self.spectral_dimensions
        ]
        sp_idx = np.arange(len(sym_path))
        indexes = [np.arange(len(item)) for item in sym_path]
        products = cartesian_product(*indexes)

        return [
            SymmetryPathway(
                channels=self.channels,
                **{
                    ch: [_ for sp, i in zip(sp_idx, item) for _ in sym_path[sp][i][ch]]
                    for ch in CHANNELS
                },
            )
            for item in products
        ]

    def _get_transition_pathways_np(self, spin_system):
        all_transitions = spin_system._all_transitions()

        isotopes = spin_system.get_isotopes(symbol=True)
        channels = [item.symbol for item in self.channels]
        if np.any([item not in isotopes for item in channels]):
            return np.asarray([])

        segments = [
            evt.filter_transitions(all_transitions, isotopes, channels)
            for dim in self.spectral_dimensions
            for evt in dim.events
            if evt.__class__.__name__ != "MixingEvent"
        ]

        segments_index = [np.arange(item.shape[0]) for item in segments]
        cartesian_index = cartesian_product(*segments_index)
        return [
            [segments[i][j] for i, j in enumerate(item)] for item in cartesian_index
        ]

    def _get_transition_pathway_weights(self, pathways, spin_system):
        if pathways == []:
            return np.asarray([])

        symbol = [item.isotope.symbol for item in spin_system.sites]
        spins = np.asarray([item.isotope.spin for item in spin_system.sites])
        channels = [item.symbol for item in self.channels]
        weights = np.ones(len(pathways), dtype=complex)

        # Mapping is a list of dict where each dict
        mapping = mixing_query_connect_map(self.spectral_dimensions)
        for obj in mapping:
            euler_angles = to_euler_list(symbol, channels, obj["mixing_query_list"])
            # Each column is list of the same angles
            alpha_, beta_, gamma_ = np.r_[euler_angles].T
            near_ = obj["near_index"]
            for j, path in enumerate(pathways):
                weights[j] *= self._calculate_transition_connect_weight(
                    path[near_[0]], path[near_[1]], spins, alpha_, beta_, gamma_
                )
        return np.round(weights, decimals=6)

    @staticmethod
    def _calculate_transition_connect_weight(trans1, trans2, spins, alpha, beta, gamma):
        """Return the transition connection weight of transitions `trans1` and `trans2`.

        Args:
            trans1: ndarray of shape (2, n_sites) representing the starting transition,
                where `n_sites` is the number of sites within the spin system. The two
                entries at axis 0 corresponds to initial (index 0) and final (index 1)
                Zeeman energy states.
            trans2: ndarray of shape (2, n_sites) representing the target transition,
                where `n_sites` is the number of sites within the spin system. The two
                entries at axis 0 corresponds to initial (index 0) and final (index 1)
                Zeeman energy states.
            spins: ndarray of spin quantum numbers of the sites.
            theta: ndarray of rotation angle per site from the mixing event.
            phi: ndarray of rotation phase per site from the mixing event.
            alpha: ndarray of first Euler angle per site from the mixing query.
            beta: ndarray of second Euler angle per site from the mixing query.
            gamma: ndarray of third Euler angle per site from the mixing query.
        """
        amp = 1
        for i, spin in enumerate(spins):
            m1_f = trans1[1][i]  # starting transition final state
            m1_i = trans1[0][i]  # starting transition initial state
            m2_f = trans2[1][i]  # landing transition final state
            m2_i = trans2[0][i]  # landing transition initial state

            amp *= transition_connect_factor(
                spin,
                m1_f,
                m1_i,
                m2_f,
                m2_i,
                alpha[i],
                beta[i],
                gamma[i],
            )
        return amp

    def _get_transition_pathway_and_weights_np(self, spin_system):
        segments = self._get_transition_pathways_np(spin_system)
        weights = self._get_transition_pathway_weights(segments, spin_system)
        return segments, weights
        # indexes = np.where(weights != 0)[0]
        # return np.asarray(segments)[indexes], weights[indexes]

    def get_transition_pathways(self, spin_system) -> List[TransitionPathway]:
        """Return a list of transition pathways from the given spin system that satisfy
        the query selection criterion of the method.

        Args:
            SpinSystem spin_system: A SpinSystem object.

        Returns:
            A list of :ref:`transition_pathway_api` objects. Each TransitionPathway
            object is an ordered collection of Transition objects.

        Example:
            >>> from mrsimulator import SpinSystem
            >>> from mrsimulator.method.lib import ThreeQ_VAS
            >>> sys = SpinSystem(sites=[{'isotope': '27Al'}, {'isotope': '29Si'}])
            >>> method = ThreeQ_VAS(channels=['27Al'])
            >>> pprint(method.get_transition_pathways(sys))
            [|1.5, -0.5⟩⟨-1.5, -0.5| ⟶ |-0.5, -0.5⟩⟨0.5, -0.5|, weight=(1+0j),
             |1.5, -0.5⟩⟨-1.5, -0.5| ⟶ |-0.5, 0.5⟩⟨0.5, 0.5|, weight=(1+0j),
             |1.5, 0.5⟩⟨-1.5, 0.5| ⟶ |-0.5, -0.5⟩⟨0.5, -0.5|, weight=(1+0j),
             |1.5, 0.5⟩⟨-1.5, 0.5| ⟶ |-0.5, 0.5⟩⟨0.5, 0.5|, weight=(1+0j)]
        """
        segments, weights = self._get_transition_pathway_and_weights_np(spin_system)
        return [
            TransitionPathway(
                [
                    Transition(initial=tr[0].tolist(), final=tr[1].tolist())
                    for tr in item
                ],
                weight=w,
            )
            for item, w in zip(segments, weights)
            if w != 0
        ]

    def _add_simple_props_to_df(self, df, prop_dict, required, drop_constant_columns):
        """Helper method for summary to reduce complexity"""
        # Iterate through property and valid Event subclass for property
        for prop, valid in prop_dict.items():
            lst = [
                getattr(ev, prop) if ev.__class__.__name__ in valid else np.nan
                for dim in self.spectral_dimensions
                for ev in dim.events
            ]
            if prop not in required and drop_constant_columns:
                # NOTE: np.isnan() cannot be passed an object (freq_contrib)
                lst_copy = np.asarray(lst)[~np.isnan(np.asarray(lst))]
                if np.unique(lst_copy).size < 2:
                    continue
            df[prop] = lst

    def summary(self, drop_constant_columns=True) -> pd.DataFrame:
        """Returns a DataFrame giving a summary of the Method. A user can specify
        optional attributes to include which appear as columns in the DataFrame. A user
        can also ask to leave out attributes which remain constant throughout the
        method. Invalid attributes for an Event will be replaced with NAN.

        Args:
            (bool) drop_constant_columns:
                Removes constant properties if True. Default is True.

        Returns:
            pd.DataFrame df:
                Event number as row and property as column. Invalid properties for an
                event type are filled with np.nan

            **Columns**

            - (str) type: Event type
            - (int) spec_dim_index: Index of spectral dimension which event belongs to
            - (str) label: Event label
            - (float) duration: Duration of the ConstantDurationEvent
            - (float) fraction: Fraction of the SpectralEvent
            - (MixingQuery) query: MixingQuery object of the MixingEvent
            - (float) magnetic_flux_density: Magnetic flux density during event in Tesla
            - (float) rotor_frequency: Rotor frequency during event in Hz
            - (float) rotor_angle: Rotor angle during event converted to Degrees
            - (FrequencyEnum) freq_contrib: Frequency

        Example:
            **All Possible Columns**

            >>> from mrsimulator.method.lib import ThreeQ_VAS
            >>> method = ThreeQ_VAS(channels=["17O"])
            >>> df = method.summary(drop_constant_columns=False)
            >>> pprint(list(df.columns))
            ['type',
             'spec_dim_index',
             'spec_dim_label',
             'label',
             'duration',
             'fraction',
             'query',
             'magnetic_flux_density',
             'rotor_frequency',
             'rotor_angle',
             'freq_contrib',
             'p',
             'd']
        """
        CD = "ConstantDurationEvent"
        SP = "SpectralEvent"
        MX = "MixingEvent"

        # Columns required to be present
        required = [
            "type",
            "label",
            "duration",
            "fraction",
            "query",
            "spec_dim_index",
            "spec_dim_label",
            "freq_contrib",
            "p",
            "d",
        ]

        # Properties which can accessed by getattr()
        prop_dict = {
            "label": (CD, SP, MX),
            "duration": CD,
            "fraction": SP,
            "query": MX,
            "magnetic_flux_density": (CD, SP),
            "rotor_frequency": (CD, SP),
            "rotor_angle": (CD, SP),
            "freq_contrib": (CD, SP),
        }

        # Create the DataFrame
        df = pd.DataFrame()

        # Populate columns which cannot be calculated from iteration
        df["type"] = [
            ev.__class__.__name__
            for dim in self.spectral_dimensions
            for ev in dim.events
        ]
        df["spec_dim_index"] = [
            i for i, dim in enumerate(self.spectral_dimensions) for _ in dim.events
        ]
        df["spec_dim_label"] = [
            dim.label for dim in self.spectral_dimensions for _ in dim.events
        ]
        self._add_simple_props_to_df(df, prop_dict, required, drop_constant_columns)

        # Add p and d symmetry pathways to dataframe
        # (future) add multi-channel support
        # IDEA: dict with "total" and "ch1"..."ch3" as keys
        df["p"] = np.transpose(
            [sym.total for sym in self.get_symmetry_pathways("P")]
        ).tolist()
        df["d"] = np.transpose(
            [sym.total for sym in self.get_symmetry_pathways("D")]
        ).tolist()

        # Convert rotor_angle to degrees
        if "rotor_angle" in df.columns:
            df["rotor_angle"] = df["rotor_angle"] * 180 / np.pi

        # Convert rotor_frequency to kHz
        if "rotor_frequency" in df.columns:
            df["rotor_frequency"] = df["rotor_frequency"] / 1000

        return df

    def plot(self, df=None, include_legend=False) -> mpl.pyplot.figure:
        """Creates a diagram representing the method. By default, only parameters which
        vary throughout the method are plotted. Figure can be finley adjusted using
        matplotlib rcParams.

        Args:
            DataFrame df:
                DataFrame to plot data from. By default DataFrame is calculated from
                summary() and will show only parameters which vary throughout the
                method plus 'p' symmetry pathway and 'd' symmetry pathway if it is not
                none or defined

            bool include_legend:
                Optional argument to include a key for event colors. Default is False
                and no key will be included in figure

        Returns:
            matplotlib.pyplot.figure

        Example:
            >>> from mrsimulator.method.lib import BlochDecaySpectrum
            >>> method = BlochDecaySpectrum(channels=["13C"])
            >>> fig = method.plot()

            **Adjusting Figure Size rcParams**

            >>> import matplotlib as mpl
            >>> from mrsimulator.method.lib import FiveQ_VAS
            >>> mpl.rcParams["figure.figsize"] = [14, 10]
            >>> mpl.rcParams["font.size"] = 14
            >>> method = FiveQ_VAS(channels=["27Al"])
            >>> fig = method.plot(include_legend=True)

            **Plotting all Parameters, including Constant**

            >>> from mrsimulator.method.lib import FiveQ_VAS
            >>> method = FiveQ_VAS(channels=["27Al"])
            >>> df = method.summary(drop_constant_columns=False)
            >>> fig = method.plot(df=df)
        """
        fig = mpl.pyplot.gcf()

        if df is None:
            df = self.summary()

        _plot(fig, df, [x.symbol for x in self.channels], include_legend)

        fig.suptitle(t=self.name if self.name is not None else "")
        return fig

    def shape(self) -> tuple:
        """The shape of the method's spectral dimension array.

        Returns:
            tuple

        Example:
            >>> from mrsimulator.method import Method
            >>> method = Method(
            ...     channels=['1H'],
            ...     spectral_dimensions=[{'count': 40}, {'count': 10}]
            ... )
            >>> method.shape()
            (40, 10)
        """
        return tuple(item.count for item in self.spectral_dimensions)
