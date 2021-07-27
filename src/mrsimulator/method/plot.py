# -*- coding: utf-8 -*-
from itertools import groupby

import matplotlib.projections as proj
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from matplotlib.patches import Rectangle
from matplotlib.ticker import FixedLocator
from matplotlib.ticker import MaxNLocator

__author__ = "Matthew D. Giammar"
__email__ = "giammar.7@buckeyemail.osu.edu"


# TODO: add globally defined alpha for rectangles
# NOTE: Matplotlib should automatically generate new colors when none specified
DURATION_WIDTH = 0.5  # Width of one ConstantDurationEvent
SPECTRAL_MULTIPLIER = 0.8  # Width multiplier for all SpectralEvents
MIXING_WIDTH = 0.25  # tip_angle of 360 degrees
COLORS = {
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
DEFAULT_ANNO_KWARGS = {
    "annotation_clip": False,
    "color": "black",
    "ha": "center",
    "va": "center",
}


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

    mix_ev:
        True where event is MixingEvent, false otherwise

    """

    name = "custom_axes"
    x_data = None
    y_data = None
    col_name = None
    xmax = None
    mix_ev = None

    def make_plot(self, x_data, y_data, col_name, mix_ev, format_kwargs, plot_kwargs):
        """Main workflow function to format and plot data on Axes"""
        self.x_data = np.array(x_data, dtype=float)
        self.y_data = np.array(y_data, dtype=float)
        self.col_name = col_name
        self.xmax = max(x_data)
        self.mix_ev = mix_ev

        self._format(xmax=self.xmax, **format_kwargs)
        self.plot(x=self.x_data, y=self.y_data, **plot_kwargs)

    def plot(self, x, y, _format=True, labels=True, **kwargs):
        """Plot formatted data"""
        if _format:
            # matplotlib automatically blanks out nan values. y may contain nans
            # mask only differs from y when plotting "rotor_frequnecy"
            x, y, y_mask = self._format_x_and_y_data(x=x, y=y)
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
                    rect_kwargs=dict(color=COLORS["inf_speed"]),
                )
        # Locate places with undefined parameters
        reigons = self._locate_reigons_with_val(x=x, y=y)
        for reigon in reigons:
            self._add_rect_with_label(
                x0=reigon[0],
                x1=reigon[1],
                label="undef",
                rect_kwargs=dict(color=COLORS["undef"]),
            )

    def _locate_reigons_with_val(self, x, y, val=np.nan):
        """Locates reigons and returns a list of tuples denoting range with value"""
        # Locate value in y
        if val is np.nan:
            loc = np.isnan(y)
        else:
            loc = y >= val

        # Find indexes where loc_k and loc_k+1 are both true
        left_x = np.argwhere(np.logical_and(loc[1:], loc[:-1])).flatten()

        # Create tuples of left and right x points
        return [(x[i], x[i + 1]) for i in left_x]

    def _format_x_and_y_data(self, x, y):
        """Removes invalid event data (nans from MixingEvents) while keeping valid but
        undefined event data in y. Extends y double to accheive stair-step pattern
        """
        y_data = y[~self.mix_ev]  # Keep data from events that are not MixingEvents
        x_data = x
        y_data = np.asarray([n for num in y_data for n in (num, num)])  # Extend double

        # Mask pseudo-infinite spinning speed with nans for plotting "rotor_frequency"
        if self.col_name == "rotor_frequency":
            mask = np.where(y_data >= 1e9, np.nan, y_data)
            return x_data, y_data, mask

        # Insert zero or remove first x point for p and d
        if self.col_name in ["p", "d"]:
            if self.mix_ev[0]:
                y_data = np.insert(y_data, 0, 0)
            else:
                x_data = x[1:]

        return x_data, y_data, y_data

    def _format(
        self,
        locator=None,
        grid={"axis": "y", "color": "black", "alpha": 0.2},
        ymargin=0.35,
        show_xaxis=False,
        show_yaxis=True,
        **kwargs,
    ):
        """Format Axes helper function

        Args:
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

        if locator is None:
            locator = MaxNLocator(nbins=4, integer=True, min_n_ticks=3)

        if self.col_name == "rotor_frequency":
            locator.set_params(integer=False, steps=[1.5, 5, 7])

        self.get_xaxis().set_visible(show_xaxis)
        self.set_xlim(0, self.xmax)

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
            y_mid = (top + bottom) / 2
            self.annotate(text=label, xy=((x1 + x0) / 2, y_mid), **anno_kwargs)


class MultiLineAxes(CustomAxes):
    """Axes subclass for multiline plots such as 'p' or 'd'"""

    name = "multi_line_axes"

    def make_plot(self, x_data, y_data, col_name, mix_ev, format_kwargs, plot_kwargs):
        """Main workflow function to format and plot data on Axes"""
        self.x_data = x_data
        self.y_data = np.stack(y_data.values).transpose()  # each row is symm pathway

        if np.asarray(self.y_data).ndim != 2:
            raise ValueError("Symmetry pathway data is misshapen. Data must be 2d")

        self.col_name = col_name
        self.xmax = max(x_data)
        self.mix_ev = mix_ev

        self._format(**format_kwargs)
        self.plot(x=self.x_data, y=self.y_data, **plot_kwargs)

    def plot(self, x, y, do_offset=True, **kwargs):
        if do_offset:
            y = self._offset_overlaps(y)

        for row in y:
            super().plot(x=x, y=row, labels=False, **kwargs)

    def _format(self, **kwargs):
        """Determines locator for axes and calls super()._format()"""
        locator = None

        # Check if y data is constant throughout events, ignoring nan
        _min = np.nanmin(self.y_data.flatten())
        if _min == np.nanmax(self.y_data.flatten()):
            self.set_ylim(_min - 1.5, _min + 1.5)
            locator = FixedLocator([_min - 1, _min, _min + 1])

        super()._format(locator=locator, **kwargs)

    def _offset_overlaps(self, y, offset_pct=0.03):
        """Offsets y at overlapping values"""
        # calculate raw offset based on percent of data range [max(y) - min(y)]
        offset = np.nanmax(self.y_data.flatten()) - np.nanmin(self.y_data.flatten())
        if offset == 0:
            offset = 2
        offset *= offset_pct

        # Check for only one symmetry pathway present
        if y.shape[0] == 1:
            return y
        data = y.transpose()  # Each row represents all symmetry pathways for an event

        # Loop over rows and find where pathways overlap
        for row in data:
            unique = np.unique(row)

            # No symmetry pathway overlap
            if unique.size == row.size:
                continue

            # Find indicies of overlap for each number and add offsets
            for num in unique:
                where = np.where(row == num)[0]
                np.add.at(
                    row,
                    where,
                    np.linspace(
                        -where.size * offset / 2,
                        where.size * offset / 2,
                        num=where.size,
                    ),
                )

        return data.transpose()


class SequenceDiagram(CustomAxes):
    """Axes subclass holding the sequence diagram of events"""

    name = "sequence_axes"
    ylim = None

    def make_plot(self, df, x_data, ylim=[0, 1]):
        """Main workflow function to plot sequence diagram"""
        self.x_data = x_data
        self.ylim = ylim
        self.xmax = max(x_data)

        self._format(grid={}, ylim=self.ylim, margin=None, axis="off")
        self._plot_sequence_diagram(df=df)

    def _format(self, **kwargs):
        # Turn of x and y axis and call super(),_format
        self.axis("off")
        self.set_ylim(*self.ylim)

        super()._format(**kwargs)

    def _format_mix_label(self, ta, p):
        """Helper method to format label for mixing events. Returns (str, float)"""
        # Format to one decimal place if float, otherwise no decimal places
        tip_angle = "{1:0.{0}f}".format(int(not float(ta).is_integer()), ta)
        phase = "{1:0.{0}f}".format(int(not float(p).is_integer()), p)
        return (f"({tip_angle}, {phase})", ta / 360 * MIXING_WIDTH)

    def _plot_spec_dims(self, df, ev_groups):
        """Adds lines and labels denoting spectral dimensions"""
        x_data = self.x_data
        ylim = self.ylim
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

    def _plot_sequence_diagram(self, df):
        """Plots sequence diagram (order of events)"""
        ev_groups = [(_type, sum(1 for _ in gp)) for _type, gp in groupby(df["type"])]

        # Last event is MixingEvent
        if ev_groups[-1][0] == "MixingEvent":
            self.x_data = np.append(self.x_data, self.x_data[-1])
            n_end_mix = ev_groups[-1][1]
            # Total angle / 360 * MIXING_WIDTH
            # Get last 'n_end_mix' events
            offset = sum(df["tip_angle"][-n_end_mix:]) / 360 * MIXING_WIDTH
            self.x_data[-2] -= offset

        self._plot_spec_dims(df=df, ev_groups=ev_groups)

        df_idx = 0
        x_idx = 0
        for _type, num in ev_groups:
            if _type == "MixingEvent":
                # Leftmost x point of next MixingEvent rectangle
                x0 = self.x_data[x_idx]
                # Iterate over each MixingEvent in group and plot rectangle
                for j in range(num):
                    label, width = self._format_mix_label(
                        ta=df["tip_angle"][df_idx + j],
                        p=df["phase"][df_idx + j],
                    )
                    super()._add_rect_with_label(
                        x0=x0,
                        x1=x0 + width,
                        label=label,
                        rect_kwargs=dict(
                            color=COLORS["MixingEvent"],
                            alpha=0.2,
                        ),
                        anno_kwargs=dict(
                            color="black",
                            ha="center",
                            va="center",
                            rotation=90,
                        ),
                    )
                    x0 += width
                x_idx += 1
            else:
                for j in range(num):
                    # Increment x_idx if no mixing event between events
                    if self.x_data[x_idx] == self.x_data[x_idx + 1]:
                        x_idx += 1
                    super()._add_rect_with_label(
                        x0=self.x_data[x_idx],
                        x1=self.x_data[x_idx + 1],
                        label=df["label"][df_idx + j],
                        rect_kwargs=dict(
                            color=COLORS[df["type"][df_idx]],
                            alpha=0.2,
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
    ]

    # Required columns are not present in df
    if not (set(required)).issubset(set(df.columns)):
        raise ValueError("Some required columns were not present in the DataFrame.")

    # Remove 'freq_contrib' since it cannot be plotted
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

    # Check if d symmetry pathway is all nan, removing if true
    if np.all(np.isnan(df["d"].tolist())):
        params.remove("d")

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


def _add_legend(fig):
    """Adds legend for event color and type"""
    fig.legend(
        handles=[
            Patch(facecolor=COLORS["MixingEvent"], alpha=0.2, label="MixingEvent"),
            Patch(facecolor=COLORS["SpectralEvent"], alpha=0.2, label="SpectralEvent"),
            Patch(
                facecolor=COLORS["ConstantDurationEvent"],
                alpha=0.2,
                label="ConstantDurationEvent",
            ),
        ],
        loc="upper left",
        fontsize="small",
    )


def _calculate_n_channels(df):
    """Calculates the number of channels present in the method DataFrame"""
    # TODO: (future) implement functionality for calcuation
    # Currently hardcoded to 1
    # Maybe move this into _format_df?
    return 1


def plot(fig, df, include_legend) -> plt.figure:
    """Plot symmetry pathways and other requested parameters on figure"""
    # TODO: (future) add functionality for multiple channels
    mix_ev = np.array(df["type"] == "MixingEvent")
    params = _format_df(df)
    x_data, x_offset = _make_normal_and_offset_x_data(df)

    # Calculations and setup gridspec object for number of subplots
    n_channels = _calculate_n_channels(df)
    nrows = n_channels + 1 + len(params)  # channels + p symmetry pathway + params
    gs = fig.add_gridspec(nrows=nrows, ncols=1)  # nrows for multiple channels
    gs_row_idx = 0

    # Register custom matplotlib projections
    proj.register_projection(CustomAxes)
    proj.register_projection(MultiLineAxes)
    proj.register_projection(SequenceDiagram)

    # Sequence diagram Axes
    seq_ax = fig.add_subplot(gs[gs_row_idx, 0], projection="sequence_axes")
    seq_ax.make_plot(df=df, x_data=x_offset)
    gs_row_idx += 1

    # p and d Axes
    p_ax = fig.add_subplot(gs[gs_row_idx, 0], projection="multi_line_axes")
    p_ax.make_plot(
        x_data=x_offset,
        y_data=df["p"],
        col_name="p",
        mix_ev=mix_ev,
        format_kwargs={},
        plot_kwargs={"alpha": 0.7},
    )
    gs_row_idx += 1

    if "d" in params:
        params.remove("d")
        d_ax = fig.add_subplot(gs[gs_row_idx, 0], projection="multi_line_axes")
        d_ax.make_plot(
            x_data=x_offset,
            y_data=df["d"],
            col_name="d",
            mix_ev=mix_ev,
            format_kwargs={},
            plot_kwargs={"alpha": 0.7},
        )
        gs_row_idx += 1

    # params Axes
    for i, param in enumerate(params, gs_row_idx):
        ax = fig.add_subplot(gs[i, 0], projection="custom_axes")
        ax.make_plot(
            x_data=x_data,
            y_data=df[param],
            col_name=param,
            mix_ev=mix_ev,
            format_kwargs={},
            plot_kwargs={},
        )

    # Add legend for event colors
    if include_legend:
        _add_legend(fig)

    return fig
