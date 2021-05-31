# -*- coding: utf-8 -*-
from copy import deepcopy
from typing import ClassVar
from typing import Dict
from typing import List
from typing import Union

import csdmpy as cp
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from mrsimulator.spin_system.isotope import Isotope
from mrsimulator.transition import Transition
from mrsimulator.transition import TransitionPathway
from mrsimulator.utils.parseable import Parseable
from pydantic import Field
from pydantic import PrivateAttr
from pydantic import validator

from .spectral_dimension import SpectralDimension
from .utils import cartesian_product

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


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

        >>> bloch = Method(channels=['1H'])
        >>> bloch.channels = ['1H']

    spectral_dimensions:
        The number of spectral dimensions depends on the given method. For example, a
        `BlochDecaySpectrum` method is a one-dimensional method and thus requires a
        single spectral dimension. The default is a single default
        :ref:`spectral_dim_api` object.

        Example
        -------

        >>> bloch = Method(channels=['1H'])
        >>> bloch.spectral_dimensions = [SpectralDimension(count=8, spectral_width=50)]
        >>> # or equivalently
        >>> bloch.spectral_dimensions = [{'count': 8, 'spectral_width': 50}]

    simulation:
        An object holding the result of the simulation. The initial value of this
        attribute is None. A value is assigned to this attribute when you run the
        simulation using the :py:meth:`~mrsimulator.Simulator.run` method.

    experiment:
        An object holding the experimental measurement for the given method, if
        available. The default value is None.

        Example
        -------

        >>> bloch.experiment = my_data # doctest: +SKIP

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

        >>> method = Method2D(channels=['1H'])
        >>> method.affine_matrix = [[1, -1], [0, 1]]
        >>> print(method.affine_matrix)
        [[1, -1], [0, 1]]
    """
    channels: List[str]
    spectral_dimensions: List[SpectralDimension] = [SpectralDimension()]
    affine_matrix: List = None
    simulation: Union[cp.CSDM, np.ndarray] = None
    experiment: Union[cp.CSDM, np.ndarray] = None

    # global
    magnetic_flux_density: float = Field(default=9.4, ge=0.0)
    rotor_frequency: float = Field(default=0.0, ge=0.0)
    rotor_angle: float = Field(default=0.9553166181245, ge=0.0, le=1.5707963268)

    _named_method: bool = PrivateAttr(False)
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
        validate_assignment = True
        arbitrary_types_allowed = True

    @validator("channels", always=True)
    def validate_channels(cls, v, *, values, **kwargs):
        return [Isotope(symbol=_) for _ in v]

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        _ = [
            setattr(ev, item, getattr(self, item))
            for sd in self.spectral_dimensions
            for ev in sd.events
            for item in self.property_units.keys()
            if hasattr(ev, item) and getattr(ev, item) is None
        ]

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
            return
        v1 = np.asarray(v)
        dim_len = len(values["spectral_dimensions"])
        if v1.size != dim_len ** 2:
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
        # mth = super().json(units=unit)
        mth = {_: self.__getattribute__(_) for _ in ["name", "label", "description"]}
        mth["channels"] = [item.json() for item in self.channels]
        mth["spectral_dimensions"] = [
            item.json(units=units) for item in self.spectral_dimensions
        ]

        # add global parameters
        evt_d = self.property_units.items()
        global_ = (
            {k: f"{self.__getattribute__(k)} {u}" for k, u in evt_d}
            if units
            else {k: self.__getattribute__(k) for k, u in evt_d}
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

    def _get_transition_pathways_np(self, spin_system):
        all_transitions = spin_system._all_transitions()

        isotopes = spin_system.get_isotopes(symbol=True)
        channels = [item.symbol for item in self.channels]
        if np.any([item not in isotopes for item in channels]):
            return []

        segments = [
            evt.filter_transitions(all_transitions, isotopes, channels)
            for dim in self.spectral_dimensions
            for evt in dim.events
        ]

        # if segments == []:
        #     return []

        segments_index = [np.arange(item.shape[0]) for item in segments]
        cartesian_index = cartesian_product(*segments_index)
        return [
            [segments[i][j] for i, j in enumerate(item)] for item in cartesian_index
        ]

    def get_transition_pathways(self, spin_system) -> List[TransitionPathway]:
        """Return a list of transition pathways from the given spin system that satisfy
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
            >>> method = Method2D(
            ...     channels=['1H'],
            ...     spectral_dimensions=[{'count': 40}, {'count': 10}]
            ... )
            >>> method.shape()
            (40, 10)
        """
        return tuple([item.count for item in self.spectral_dimensions])

    def events_to_dataframe(self, properties) -> pd.DataFrame:
        """Returns dataframe of requested Event properties

        Args:
            List properties: list of properties to keep

        Returns:
            pd.DataFrame df: properties as columns and event number as row
        """
        data = [e.json(units=False) for d in self.spectral_dimensions for e in d.events]
        df = pd.DataFrame(data)
        # print(df)
        # NOTE: Should p and d be capatalized
        P_list = np.array(
            [
                [
                    sum([sum(q[ch]["P"]) if ch in q else 0 for q in tq])
                    for ch in ["ch1", "ch2", "ch3"]
                ]
                for tq in df["transition_query"]
            ]
        )
        D_list = np.array(
            [
                [
                    sum([sum(q[ch]["D"]) if ch in q else 0 for q in tq])
                    for ch in ["ch1", "ch2", "ch3"]
                ]
                for tq in df["transition_query"]
            ]
        )
        # Remove channels which are all zero
        if not np.any(P_list[:, 2]) and not np.any(D_list[:, 2]):
            P_list = np.delete(P_list, 2, 1)
            D_list = np.delete(D_list, 2, 1)
        if not np.any(P_list[:, 1]) and not np.any(D_list[:, 1]):
            P_list = np.delete(P_list, 1, 1)
            D_List = np.delete(D_list, 1, 1)
        df["p"] = P_list.tolist()
        df["d"] = D_list.tolist()

        # (not properties) intersect (columns) + 'fraction' and 'duration
        keep = set(properties + ["fraction", "duration"])
        drop = set(df.columns).symmetric_difference(keep).intersection(set(df.columns))
        # df = df.where(pd.notnull(df), None)
        return df.drop(drop, axis=1)

    def _make_x_points(self, df):
        """Returns list of x points to use in plotting"""
        # TODO: Find value for constant duration event (tau/2)?
        # Question: Will method always have a SpectralEvent
        points = [0]
        dur_present = "duration" in df.columns


        DURATION_WIDTH = 0.4

        for i, row in df.iterrows():
            if row["fraction"] is not None:
                # SpectralEvent
                # Question: Do SpectralEvents have a total time to add up to
                next_x = points[-1] + row["fraction"]
                points.extend((next_x, next_x))

            elif dur_present and row["duration"] is not None:
                # ConstantDurationEvent
                next_x = points[-1] + DURATION_WIDTH
                points.extend((next_x, next_x))

            elif row["fraction"] is None and dur_present and row["duration"] is None:
                # MixingEvent
                # TODO: implement mixing event widths
                # - cut out time from neiboriung event
                # - account for multiple adjacent mixing events
                points.extend((points[-1], points[-1]))

        points.pop()
        return points

    def _plot_on_axis(self, ax, x_data, y_data, name):
        """Helper method to plot data and do formatting on axes"""
        if name == "p" or name == "d":
            # y_data is a list of lists by chanel number
            colors = ["red", "green", "blue"]
            for ch in range(len(y_data[0])):
                ax.plot(
                    x_data,
                    [_ for num in y_data[ch] for _ in (num, num)],
                    color=colors.pop(),
                    alpha=0.6,
                )
            # TODO: Clean up formatting of y-ticks
            _min = min([num for tup in y_data for num in tup])
            _max = max([num for tup in y_data for num in tup])
            _range = max(abs(_min), abs(_max))
            step = 1 if _range < 5 else 2
            ax.set_ylim(-(_range+step), _range+step)
            ax.set_yticks(np.arange(-_range, _range+1, step=step))

            ax.grid(axis="y", color="black", alpha=0.2)
        else:
            # y_data is a list of numbers
            ax.plot(x_data, [_ for num in y_data for _ in (num, num)], color="b", alpha=0.6)
        # y-axis formatting
        ax.set_ylabel(name.replace("_", " "))
        # Can y-axis numbers format better
        ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f"))
        for side in ["top", "bottom", "left", "right"]:
            ax.spines[side].set_linewidth(1.5)

    def plot(self, properties=["p", "d"]) -> plt.Figure:
        """Create Plotly figure of symmetry pathways for method

        Args:
            List properties: list properties to plot. Default p and d transitions

        Returns:
            Figure fig: Matplotlib figure

        Example:
            TODO add example code
        """
        # Should p and d be capatalized?

        if properties is None:
            return plt.figure()

        df = self.events_to_dataframe(properties)
        fig, axs = plt.subplots(
            nrows=len(properties) + 1,
            ncols=1,
            figsize=(df.shape[0]*2, len(properties)*1.5 + 2), 
            sharex=True,
            gridspec_kw={"hspace": 0.0},
        )
        # axs should always be longer than 1 element, otherwise cast as list

        x_points = self._make_x_points(df)
        print(x_points)
        axs[0].set_xlim(0, x_points[-1])
        # need x ticks?

        # Setup pulse sequence diagram
        def exp_sin(x):
            # NOTE: Should array be zero anchored?
            if x[0] != 0:
                x = x - x[0]
            return np.exp(-10*x) * np.sin(28*np.pi*x)
        for i, row in df.iterrows():
            dur_present = "duration" in df.columns
            if row["fraction"] is not None:
                # SpectralEvent
                arr = np.linspace(x_points[i*2], x_points[(i*2)+1], 100)
                axs[0].plot(arr, exp_sin(arr), color="black")
            elif dur_present and row["duration"] is not None:
                # ConstantDurationEvent
                axs[0].text(x=x_points[i*2], y=0.5, s="tau/2")
            elif row["fraction"] is None and dur_present and row["duration"] is None:
                # MixingEvent
                axs[0].plot([x_points[i*2], x_points[(i*2)+1]], [0, 1], color="black")
        axs[0].set_axis_off()

        # Iterate through axes and plot data
        for i, ax in enumerate(axs[1:], 0):
            self._plot_on_axis(ax, x_points, df[properties[i]], properties[i])

        # TODO: Add color key per channel

        fig.suptitle(self.name if self.name is not None else "Name")
        fig.tight_layout()
        return fig
