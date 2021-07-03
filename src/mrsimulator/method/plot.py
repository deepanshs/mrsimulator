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
DURATION_WIDTH = 0.5  # Width of one ConstantDurationEvent
SPECTRAL_MULTIPLIER = 0.8  # Width multiplier for all SpectralEvents
MIXING_WIDTH = 0.25  # tip_angle of 360 degrees
EVENT_COLORS = {
    "ConstantDurationEvent": "orange",
    "SpectralEvent": "g",
    "MixingEvent": "b",
    "inf_speed": "y",
    "nan": "k",
}
LABLES = {
    "rotor_angle": r"$\theta_r$ / deg",
    "rotor_frequency": r"$\nu_r$ / kHz",
    "magnetic_flux_density": r"$B_0$ / T",
}
DEFAULT_ANNO_KWARGS = dict(
    annotation_clip=False,
    color="black",
    ha="center",
    va="center",
    fontsize=8,
)


class CustomAxes(plt.Axes):
    """Subclass of matplotlib Axes holding spesific funtions and data

    Attributes
    ----------

    name:
        String of the matplotlib projection name of this Axes subclass. This value
        is not intended to be changed

    x_data:
        Unformatted x-axis data to be plotted

    y_data:
        Unformatted y-axis data to be plotted

    col_name:
        Column which the y_data came from

    xmax:
        Maximum x value in x_data. Used to set xlim to (0, xmax)

    """

    name = "custom_axes"
    x_data = None
    y_data = None
    col_name = None
    xmax = None

    def make_plot(self, x_data, y_data, col_name, format_kwargs, plot_kwargs):
        """Main workflow function to format and plot data on Axes"""
        self.x_data = x_data
        self.y_data = y_data
        self.col_name = col_name

        xmax = max(x_data)

        self._format(xmax=xmax, **format_kwargs)
        self.plot(x=x_data, y=y_data, **plot_kwargs)

    def plot(self, x, y, _format=True, labels=True, mix_ev=[], **kwargs):
        """Plot formatted data"""
        if _format:
            # matplotlib automatically blanks out nan values. y may contain nans
            # mask only differs from y when plotting "rotor_frequnecy"
            x, y, y_mask = self._format_x_and_y_data(x=x, y=y, mix_ev=mix_ev)
        else:
            y_mask = y

        # 'x' and 'y' cannot be passed as keyword args to plt.Axes.plot()
        super().plot(x, y_mask, **kwargs)

        # labels must come after plot to setup y-axis limits
        if labels:
            self._add_blank_space_labels(x=x, y=y)

    def _add_blank_space_labels(self, x, y):
        """Adds labels for blank spaces in plot"""
        if self.col_name == "rotor_frequency":
            # Locate places with infinite spinning speed (1e9 Hz)
            reigons = self._locate_reigons_with_val(x=x, y=y, val=1e9)
            for reigon in reigons:
                self._add_rect_with_label(
                    x0=reigon[0],
                    x1=reigon[1],
                    label="inf speed",
                    rect_kwargs=dict(color=EVENT_COLORS["inf_speed"]),
                )
        # Locate places with undefined parameters
        reigons = self._locate_reigons_with_val(x=x, y=y)
        for reigon in reigons:
            self._add_rect_with_label(
                x0=reigon[0],
                x1=reigon[1],
                label="undef",
                rect_kwargs=dict(color=EVENT_COLORS["undef"]),
            )

    def _locate_reigons_with_val(self, x, y, val=np.nan):
        """Locates reigons and returns a list of tuples denoting range with value"""
        # Locate value in y
        y = np.asarray(y)
        if val is np.nan:
            loc = np.isnan(y)
        else:
            loc = y >= val

        # Find indexes where loc_k and loc_k+1 are both true
        left_x = np.argwhere(np.logical_and(loc[1:], loc[:-1])).flatten()

        # Create tuples of left and right x points
        return [(x[i], x[i + 1]) for i in left_x]

    def _format_x_and_y_data(self, x, y, mix_ev):
        """Removes invalid event data (nans from MixingEvents) while keeping valid but
        undefined event data in y. Extends y double to accheive stair-step pattern
        """
        y_data = y[~mix_ev]  # Keep data from events that are not MixingEvents
        x_data = x
        y_data = np.asarray([n for num in y_data for n in (num, num)])  # Extend double

        # Mask pseudo-infinite spinning speed with nans for plotting "rotor_frequency"
        if self.col_name == "rotor_frequency":
            mask = np.where(y_data >= 1e9, np.nan, y_data)
            return x_data, y_data, mask

        # Insert zero or remove first x point for p and d
        if self.col_name in ["p", "d"]:
            if mix_ev[0]:
                y_data = np.insert(y_data, 0, 0)
            else:
                x_data = x[1:]

        return x_data, y_data, y_data

    def _calculate_range(self, data):
        """Helper function to calculate range of given data"""
        return (np.nanmin(data), np.nanmax(data))  # min, max ignoring nans

    def _format(
        self,
        xmax=None,
        locator=MaxNLocator(nbins=4, integer=True, min_n_ticks=3),
        grid={"axis": "y", "color": "black", "alpha": 0.2},
        ymargin=0.35,
        show_xaxis=False,
        show_yaxis=True,
        **kwargs
    ):
        """Format Axes helper function

        Args:
            (float) xmax:
                x-axis limit to set

            (matplotlib.ticker) locator:
                Object which determines the ticks on the y-axis

            (dict) grid:
                Grid kwargs which determine the appearance of the grid

            (float) ymargin:
                Value to be passed to matplotlib.axes.Axes.margin for padding

            (bool) show_xaxis:
                Shows x-axis and ticks if True, otherwise hidden

            (bool) show_yaxis:
                Shows y-axis and ticks if True, otherwise hidden
        """
        label = LABLES[self.col_name] if self.col_name in LABLES else self.col_name
        if self.col_name == "rotor_frequency":
            locator.set_params(integer=False, steps=[1.5, 5, 7], nbins=5)
            pass

        self.get_xaxis().set_visible(show_xaxis)
        self.set_xlim(0, xmax)

        self.get_yaxis().set_visible(show_yaxis)
        self.get_yaxis().set_major_locator(locator)
        self.set_ylabel(label)

        self.grid(**grid)
        self.margins(y=ymargin)

    def _add_rect_with_label(self, x0, x1, label, rect_kwargs, anno_kwargs={}):
        """Add a rectangle between x0 and x1 on ax representing event"""
        # [height, color, alpha] required in rect_kwargs
        bottom, top = self.get_ylim()
        if "height" not in rect_kwargs:
            rect_kwargs["height"] = top - bottom

        if "alpha" not in rect_kwargs:
            rect_kwargs["alpha"] = 0.2

        if "color" not in rect_kwargs:
            raise ValueError("No color in `rect_kwargs`. A color must be spesified")

        if anno_kwargs == {}:
            anno_kwargs = DEFAULT_ANNO_KWARGS

        rect = Rectangle(xy=(x0, bottom), width=x1 - x0, **rect_kwargs)
        self.add_patch(rect)
        if label is not None:
            y_mid = (top - bottom) / 2
            self.annotate(text=label, xy=((x1 + x0 + 0.025) / 2, y_mid), **anno_kwargs)


class MultiLineAxes(CustomAxes):
    """Axes subclass for multiline plots such as 'p' or 'd'"""

    name = "multi_line_axes"

    def make_plot(self, x_data, y_data, col_name, format_kwargs, plot_kwargs):
        """Main workflow function to format and plot data on Axes"""
        self.x_data = x_data
        self.y_data = y_data
        self.col_name = col_name

        xmax = max(x_data)
        # Cast y_data to numpy array (should be 2d)
        y_data = np.stack(y_data.values)
        if np.asarray(y_data).ndim != 2:
            raise ValueError("Symmetry pathway data is misshapen. Data must be 2d")

        self._format(
            xmax=xmax,
            locator=MaxNLocator(nbins=5, steps=[1, 2, 3], integer=True, min_n_ticks=3),
            **format_kwargs,
        )
        self.plot(x=x_data, y=y_data, **plot_kwargs)

    def plot(self, x, y, **kwargs):
        for i in range(len(y[0])):
            # ith column (ith symmetry pathway)
            super().plot(x=x, y=y[:, i], labels=False, **kwargs)


class SequenceDiagram(CustomAxes):
    """Axes subclass holding the sequence diagram of events"""

    name = "sequence_axes"

    # TODO: Fix formatting of spectral dimension labels
    # TODO: Implement dynamic height based on mixing event labels

    def make_plot(self, df, x_data):
        """Main workflow function to plot sequence diagram"""
        xmax = max(x_data)
        ylim = [0, 1]  # Adjust y limits of plot here

        self._format(xmax=xmax, grid={}, ylim=ylim, margin=None, axis="off")
        self._plot_sequence_diagram(df=df, x_data=x_data, ylim=ylim)

    def _format(self, ylim=None, **kwargs):
        # Turn of x and y axis and call super(),_format
        self.axis("off")

        if ylim is not None:
            self.set_ylim(*ylim)

        super()._format(**kwargs)

    def _format_mix_label(self, tip_angle, phase):
        """Helper method to format label for mixing events. Returns (str, float)"""
        return (
            "({0:.1f}, {1:.1f})".format(tip_angle, phase),
            tip_angle / 360 * MIXING_WIDTH,
        )

    def _plot_spec_dims(self, df, x_data, ylim, ev_groups):
        """Adds lines and labels denoting spectral dimensions"""
        plot_kwargs = dict(_format=False, labels=False)

        self.plot(x=[0, 0], y=ylim, color="black", **plot_kwargs)

        last_spec_dim_x = 0
        spec_dim_idx = 0
        df_idx = 0
        x_idx = 0
        for _type, num in ev_groups:
            for j in range(num):
                if x_idx < len(x_data) - 1 and x_data[x_idx] == x_data[x_idx + 1]:
                    x_idx += 1
                if df["spec_dim_index"][df_idx + j] != spec_dim_idx:  # Next spec dim
                    x0 = x_data[x_idx]
                    if df["type"][df_idx + j] == "MixingEvent":
                        x0 += sum(df["tip_angle"][df_idx:j]) / 360.0 * MIXING_WIDTH
                    self.plot(x=[x0, x0], y=ylim, color="black", **plot_kwargs)
                    self.annotate(
                        text=df["spec_dim_label"][df_idx + j - 1],
                        xy=((last_spec_dim_x + x0) / 2, ylim[1] + 0.1),
                        **DEFAULT_ANNO_KWARGS,
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
            **DEFAULT_ANNO_KWARGS,
        )

    def _plot_sequence_diagram(self, df, x_data, ylim):
        """Plots sequence diagram (order of events)"""
        ev_groups = [(_type, sum(1 for _ in gp)) for _type, gp in groupby(df["type"])]

        if ev_groups[-1][0] == "MixingEvent":
            x_data = np.append(x_data, x_data[-1])
            # Total angle / 360 * MIXING_WIDTH
            gp__ = -ev_groups[0][1]
            offset = sum(df["tip_angle"][gp__:]) / 180 * MIXING_WIDTH
            x_data[-2] -= offset

        self._plot_spec_dims(df=df, x_data=x_data, ylim=ylim, ev_groups=ev_groups)

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
                        rect_kwargs=dict(
                            color=EVENT_COLORS["MixingEvent"],
                            alpha=0.2,
                            # height=ylim[1] - ylim[0],
                        ),
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
                            color=EVENT_COLORS[df["type"][df_idx]],
                            alpha=0.2,
                            # height=ylim[1] - ylim[0],
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
    # TODO: add kwargs [include_key, figsize, hspace] for user control
    # TODO: add scale size for plot scalability
    # TODO: (future) add functionality for multiple channels
    mix_ev = np.array(df["type"] == "MixingEvent")
    params = _format_df(df)
    x_data, x_offset = _make_normal_and_offset_x_data(df)

    # Register custom matplotlib projections
    proj.register_projection(CustomAxes)
    proj.register_projection(MultiLineAxes)
    proj.register_projection(SequenceDiagram)

    fig = plt.figure(figsize=[10, 7.5])  # Adjust figure size here
    gs = fig.add_gridspec(nrows=len(params) + 3, ncols=1, hspace=0.25)

    # Sequence diagram Axes
    seq_ax = fig.add_subplot(gs[0, 0], projection="sequence_axes")
    seq_ax.make_plot(df=df, x_data=x_offset)

    # p and d Axes
    p_ax = fig.add_subplot(gs[1, 0], projection="multi_line_axes")
    p_ax.make_plot(
        x_data=x_offset,
        y_data=df["p"],
        col_name="p",
        format_kwargs={},
        plot_kwargs=dict(mix_ev=mix_ev),
    )
    d_ax = fig.add_subplot(gs[2, 0], projection="multi_line_axes")
    d_ax.make_plot(
        x_data=x_offset,
        y_data=df["d"],
        col_name="d",
        format_kwargs={},
        plot_kwargs=dict(mix_ev=mix_ev),
    )

    # params Axes
    for i, param in enumerate(params, 3):
        ax = fig.add_subplot(gs[i, 0], projection="custom_axes")
        ax.make_plot(
            x_data=x_data,
            y_data=df[param],
            col_name=param,
            format_kwargs={},
            plot_kwargs=dict(mix_ev=mix_ev),
        )

    return fig
