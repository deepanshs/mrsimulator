# -*- coding: utf-8 -*-
from itertools import groupby

import matplotlib.projections as proj
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator

__author__ = "Matthew D. Giammar"
__email__ = "giammar.7@buckeyemail.osu.edu"


# NOTE: Matplotlib should automatically generate new colors when none specified
# NOTE: Cannot force only even ticks with MaxNLocator
DURATION_WIDTH = 0.6  # Width of one ConstantDurationEvent
SPECTRAL_MULTIPLIER = 1  # Width multiplier for all SpectralEvents
MIXING_WIDTH = 0.3  # tip_angle of 360 degrees
EVENT_COLORS = {  # TODO: add colors
    "ConstantDurationEvent": "orange",
    "SpectralEvent": "g",
    "MixingEvent": "b",
}
LABLES = {
    "rotor_angle": r"$\theta_r$ / deg",
    "rotor_frequency": r"$\nu_r$ / kHz",
    "magnetic_flux_density": r"$B_0$ / T",
}


class CustomAxes(plt.Axes):
    """Subclass of matplotlib Axes holding spesific funtions and data"""

    name = "custom_axes"
    locator = MaxNLocator(nbins=5, integer=True, min_n_ticks=3)
    _grid = {"axis": "y", "color": "black", "alpha": 0.2}

    def make_plot(self, x_data, y_data, col_name, format_kwargs, plot_kwargs):
        """Main workflow function to format and plot data on Axes"""
        self._format(col_name=col_name, **format_kwargs)
        self.plot(x=x_data, y=y_data, **plot_kwargs)

    def plot(self, x, y, _format=True, labels=True, first_ev_mix=None, **kwargs):
        """Plot formatted data"""
        if _format:
            x, y, nans = self._format_x_and_y_data(x=x, y=y, first_ev_mix=first_ev_mix)
        else:
            nans = None

        # TODO: Make blank spaces in line

        if labels:
            # Add labels for blank spaces
            # NOTE: Check for nans is none
            # TODO: create labels with self._add_rect_with_label
            # TODO: Loop over labels
            self._add_blank_space_labels(x_data=x, nans=nans)
            pass

        # 'x' and 'y' cannot be passed as keyword args to plt.Axes.plot()
        super().plot(x, y, **kwargs)

    def _add_blank_space_labels(self, x_data, nans, col_name=None):
        """Adds labels for blank spaces in plot"""
        # Computation
        # Call self._add_rect_with_label
        pass

    def _format_x_and_y_data(self, x, y, first_ev_mix=None):
        """Formats x & y data for plotting and locates NANs"""
        x_data = x
        nans = np.isnan(y)
        y_data = y[~nans]  # Clean nans from data
        y_data = [n for num in y_data for n in (num, num)]  # Extend y_data double

        # Insert zero or remove first x point
        if first_ev_mix:
            y_data = [0] + y_data
        elif first_ev_mix is False:
            x_data = x[1:]

        return x_data, y_data, nans

    def _format(self, col_name, **kwargs):
        """Format Axes helper function"""
        label = LABLES[col_name] if col_name in LABLES else col_name
        if col_name == "rotor_angle":
            # TODO: Create new locator to clean up y_ticks
            pass

        if "locator" in kwargs:
            self.locator = kwargs["locator"]

        if "_grid" in kwargs:
            self._grid = kwargs["_grid"]

        if "y_lim" in kwargs:
            self.set_ylim(*(kwargs["ylim"]))

        if "axis" in kwargs:
            self.axis(kwargs["axis"])

        self.get_xaxis().set_visible(False)
        # set x_lims
        self.set_ylabel(label)
        self.get_yaxis().set_major_locator(self.locator)
        # setup y_limits with space buffer
        self.grid(**self._grid)

    def _add_rect_with_label(self, x0, x1, label, rect_kwargs, anno_kwargs):
        """Add a rectangle between x0 and x1 on ax representing event"""
        # [height, color, alpha] required in rect_kwargs
        # [color, ha, va, fontsize, rotation] required in anno_kwargs
        if "height" not in rect_kwargs:
            bottom, top = self.get_ylim()
            rect_kwargs["height"] = top - bottom

        rect = Rectangle(xy=(x0, 0), width=x1 - x0, **rect_kwargs)
        self.add_patch(rect)
        if label is not None:
            self.annotate(text=label, xy=((x1 + x0 + 0.025) / 2, 0.5), **anno_kwargs)


class MultiLineAxes(CustomAxes):
    """Axes subclass for multiline plots such as 'p' or 'd'"""

    name = "multi_line_axes"
    locator = MaxNLocator(nbins=5, steps=[1, 2, 3], integer=True, min_n_ticks=3)

    def make_plot(self, x_data, y_data, col_name, format_kwargs, plot_kwargs):
        """Main workflow function to format and plot data on Axes"""
        # Cast y_data to numpy array (should be 2d)
        y_data = np.stack(y_data.values)
        if np.asarray(y_data).ndim != 2:
            raise ValueError("Symmetry pathway data is misshapen. Data must be 2d")

        super()._format(col_name=col_name, locator=self.locator, **format_kwargs)
        self.plot(x=x_data, y=y_data, **plot_kwargs)

    def plot(self, x, y, **kwargs):
        # Loop through symmetry pathways and plot
        for i in range(len(y[0])):
            # ith column (ith symmetry pathway)
            y_data = [0] + y

            super().plot(x=x, y=y_data, labels=False, **kwargs)


class SequenceDiagram(CustomAxes):
    """Axes subclass holding the sequence diagram of events"""

    name = "sequence_axes"

    def make_plot(self, df, x_data):
        """Main workflow function to plot sequence diagram"""
        super()._format(col_name="", **dict(_grid={}, axis="off"))
        self._plot_sequence_diagram(df=df, x_data=x_data)

    def _format_mix_label(self, tip_angle, phase):
        """Helper method to format label for mixing events. Returns (str, float)"""
        return (
            "({0:.1f}, {1:.1f})".format(tip_angle, phase),
            tip_angle / 360 * MIXING_WIDTH,
        )

    def _plot_spec_dims(self, df, x_data, ev_groups):
        """Adds lines and labels denoting spectral dimensions"""
        plot_kwargs = dict(_format=False, labels=False)
        anno_kwargs = dict(
            annotation_clip=False,
            color="black",
            ha="center",
            va="center",
            fontsize=8,
        )
        ylim = self.get_ylim()

        self.plot(x=[0, 0], y=ylim, color="black", **plot_kwargs)

        last_spec_dim_x = 0
        spec_dim_idx = 0
        df_idx = 0
        x_idx = 0
        for _type, num in ev_groups:
            for j in range(num):
                if x_data[x_idx] == x_data[x_idx + 1]:
                    x_idx += 1
                if df["spec_dim_index"][df_idx + j] != spec_dim_idx:  # Next spec dim
                    x0 = x_data[x_idx]
                    if df["type"][df_idx + j] == "MixingEvent":
                        x0 += sum(df["tip_angle"][df_idx:j]) / 360.0 * MIXING_WIDTH
                    self.plot(x=[x0, x0], y=ylim, color="black", **plot_kwargs)
                    self.annotate(
                        text=df["spec_dim_label"][df_idx + j - 1],
                        xy=((last_spec_dim_x + x0) / 2, ylim[1] + 0.1),
                        **anno_kwargs,
                    )
                    last_spec_dim_x = x0
                    spec_dim_idx += 1
                x_idx += 1
            df_idx += num

        # Plot last spectral dimension
        self.plot(x=[max(x_data)] * 2, y=ylim, color="black", **plot_kwargs)
        self.annotate(
            text=df["spec_dim_label"][df_idx - 1],
            xy=((last_spec_dim_x + max(x_data)) / 2, ylim[1] + 0.1),
            **anno_kwargs,
        )

    def _plot_sequence_diagram(self, df, x_data):
        """Plots sequence diagram (order of events)"""
        ev_groups = [(_type, sum(1 for _ in gp)) for _type, gp in groupby(df["type"])]

        if ev_groups[-1][0] == "MixingEvent":
            x_data = np.append(x_data, x_data[-1])
            # Total angle / 360 * MIXING_WIDTH
            gp__ = -ev_groups[0][1]
            offset = sum(df["tip_angle"][gp__:]) / 180 * MIXING_WIDTH
            x_data[-2] -= offset

        self._plot_spec_dims(df=df, x_data=x_data, ev_groups=ev_groups)

        df_idx = 0
        x_idx = 0
        for _type, num in ev_groups:
            if _type == "MixingEvent":
                # Leftmost x point of next MixingEvent rectangle
                x0 = x_data[x_idx]
                # Iterate over each MixingEvent in group and plot rectangle
                for j in range(num):
                    label, width = self._format_mix_label(
                        tip_angle=df["tip_angle"][df_idx + j],
                        phase=df["phase"][df_idx + j],
                    )
                    super()._add_rect_with_label(
                        x0=x0,
                        x1=x0 + width,
                        label=label,
                        rect_kwargs=dict(color=EVENT_COLORS["MixingEvent"], alpha=0.2),
                        anno_kwargs=dict(
                            color="black",
                            ha="center",
                            va="center",
                            fontsize=7,
                            rotation=90,
                        ),
                    )
                    x0 += width
                x_idx += 1
            else:
                for j in range(num):
                    # Increment x_idx if no mixing event between events
                    if x_data[x_idx] == x_data[x_idx + 1]:
                        x_idx += 1
                    super()._add_rect_with_label(
                        x0=x_data[x_idx],
                        x1=x_data[x_idx + 1],
                        label=df["label"][df_idx + j],
                        rect_kwargs=dict(
                            color=EVENT_COLORS[df["type"][df_idx]], alpha=0.2
                        ),
                        anno_kwargs=dict(
                            color="black",
                            ha="center",
                            va="center",
                            fontsize=7,
                        ),
                    )
                    x_idx += 1

            df_idx += num


def _check_columns(df):
    """Helper method to ensure required columns are present. Returns list of included
    optional parameters in dataframe columns.
    """
    required = [
        "type",
        "label",
        "spec_dim_index",
        "spec_dim_label",
        "duration",
        "fraction",
        "mixing_query",
        "p",
        "d",
    ]

    # Required columns are not present in df
    if not (set(required)).issubset(set(df.columns)):
        raise ValueError("Some required columns were not present in the DataFrame.")

    if "freq_contrib" in df.columns:
        required += ["freq_contrib"]

    return list(df.columns.drop(required))


def _add_tip_angle_and_phase(df):
    """Add tip_angle and phase columns to dataframe from mixing_query"""
    # NOTE Only columns for ch1 are created (1 ch only for now)
    # NOTE What should empty MixingQuerys add? (MixingQuery default?)
    # BUG: Empty mixing queries (i.e. "mixing_query": {}) gives query of `None` causing
    #       error to be thrown
    df["tip_angle"] = [
        query.ch1.tip_angle * 180 / np.pi
        if query.__class__.__name__ == "MixingQuery"
        else np.nan
        for query in df["mixing_query"]
    ]
    df["phase"] = [
        query.ch1.phase * 180 / np.pi
        if query.__class__.__name__ == "MixingQuery"
        else np.nan
        for query in df["mixing_query"]
    ]


def _format_df(df):
    """Formats dataframe and returns params to plot"""
    params = _check_columns(df)
    _add_tip_angle_and_phase(df)
    return params


def _make_x_data(df):
    """Returns list of x points to use in plotting"""
    points = [0]

    for i, row in df.iterrows():
        if row["type"] == "SpectralEvent":
            next_x = points[-1] + (row["fraction"] * SPECTRAL_MULTIPLIER)
            points.extend((next_x, next_x))
        elif row["type"] == "ConstantDurationEvent":
            next_x = points[-1] + DURATION_WIDTH
            points.extend((next_x, next_x))

    points.pop()
    return points


def _offset_x_data(df, x_data):
    """Offsets x_data based on MixingEvents"""
    offset_x = np.array([0] + x_data)
    ev_groups = [(_type, sum(1 for _ in group)) for _type, group in groupby(df["type"])]

    # Remove MixingEvents from end of sequence
    if ev_groups[-1][0] == "MixingEvent":
        ev_groups.pop()

    df_idx = 0
    # Extend first jump if first event(s) are MixingEvent
    if ev_groups[0][0] == "MixingEvent":
        gp__ = ev_groups[0][1]
        offset = sum(df["tip_angle"][0:gp__]) / 360.0 * MIXING_WIDTH
        offset_x[1] += offset
        # Increment event indexer by number of MixingEvents in first group
        df_idx += gp__
        ev_groups.pop(0)

    x_idx = 1
    for _type, num in ev_groups:
        if _type == "MixingEvent":
            up_lim__ = df_idx + num
            offset = sum(df["tip_angle"][df_idx:up_lim__]) / 360.0 * MIXING_WIDTH
            offset_x[x_idx] -= offset / 2
            offset_x[x_idx + 1] += offset / 2
            x_idx += 1
        else:
            if offset_x[x_idx] == offset_x[x_idx + 1]:
                x_idx += 1
            x_idx += (num * 2) - 1
        # Increment event indexer by number of Events in group
        df_idx += num

    return offset_x


def _make_normal_and_offset_x_data(df):
    """Calculates proper x_data and returns normal and offset x_data"""
    x_data = _make_x_data(df)

    if len(x_data) == 0:
        raise ValueError(
            "The DataFrame does not contain any SpectralEvents or "
            "ConstandDurationEvents. At least one must be present to construct a plot"
        )

    return x_data, _offset_x_data(df, x_data)


def plot(df) -> plt.figure:
    """Create figure of symmetry pathways for DataFrame representation of method"""
    # TODO: (future) add functionality for multiple channels
    first_ev_mix = df["type"][0] == "MixingEvent"
    params = _format_df(df)
    x_data, x_offset = _make_normal_and_offset_x_data(df)

    # Register custom matplotlib projections
    proj.register_projection(CustomAxes)
    proj.register_projection(MultiLineAxes)
    proj.register_projection(SequenceDiagram)

    fig = plt.figure(figsize=[6.4, 4.8])  # Adjust figure size here
    gs = fig.add_gridspec(nrows=len(params) + 3, ncols=1, hspace=0.15)

    # Sequence diagram Axes
    print("sequence")
    seq_ax = fig.add_subplot(gs[0, 0], projection="sequence_axes")
    seq_ax.make_plot(df=df, x_data=x_offset)

    # p and d Axes
    print("p")
    p_ax = fig.add_subplot(gs[1, 0], projection="multi_line_axes")
    p_ax.make_plot(
        x_data=x_offset,
        y_data=df["p"],
        col_name="p",
        format_kwargs={},
        plot_kwargs=dict(first_ev_mix=first_ev_mix),
    )
    print("d")
    d_ax = fig.add_subplot(gs[2, 0], projection="multi_line_axes")
    d_ax.make_plot(
        x_data=x_offset,
        y_data=df["d"],
        col_name="d",
        format_kwargs={},
        plot_kwargs=dict(first_ev_mix=first_ev_mix),
    )

    # params Axes
    for i, param in enumerate(params, 3):
        ax = fig.add_subplot(gs[i, 0], projection="custom_axes")
        ax.make_plot(
            x_data=x_data,
            y_data=df[param],
            col_name=param,
            format_kwargs={},
            plot_kwargs={},
        )

    return fig
