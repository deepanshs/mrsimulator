# -*- coding: utf-8 -*-
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

__author__ = "Matthew D. Giammar"
__email__ = "giammar.7@buckeyemail.osu.edu"


DURATION_WIDTH = 0.6
MIXING_WIDTH = 0.04
# TODO: Ensure cannot run out of colors
COLORS = list(mcolors.TABLEAU_COLORS)
EVENT_COLORS = {  # TODO: add colors
    "ConstantDurationEvent": 'k',
    "SpectralEvent": 'g',
    "MixingEvent": 'b'
}


def _make_x_data(df):
    """Returns list of x points to use in plotting"""
    # TODO: Find value for constant duration event (tau/2)?
    points = [0]

    for i, row in df.iterrows():
        if row["type"] == "SpectralEvent":
            next_x = points[-1] + row["fraction"]
            points.extend((next_x, next_x))
        elif row["type"] == "ConstantDurationEvent":
            next_x = points[-1] + DURATION_WIDTH
            points.extend((next_x, next_x))

    points.pop()
    return points


def _offset_x_data(df, x_data):
    """Offsets x_data based on MixingEvents"""
    # NOTE: Mixing event at end of sequence will not be concidered
    offset_x = [0] + x_data
    first_non_mix_ev_idx = np.where(df["type"] != "MixingEvent")[0][0]

    # Mixing event at begining
    if df["type"][0] == "MixingEvent":
        offset_x[1] += MIXING_WIDTH * 2 * first_non_mix_ev_idx

    mix_count = 0
    idx = 0
    for ev_type in df["type"][first_non_mix_ev_idx:]:
        if ev_type == "MixingEvent":
            mix_count += 1
            continue
        if mix_count != 0:
            offset_x[idx * 2] -= MIXING_WIDTH * mix_count
            offset_x[idx * 2 + 1] += MIXING_WIDTH * mix_count
        idx += 1
        mix_count = 0

    return offset_x


def _plot_sequence_diagram(ax, x_data, df):
    """Helper method to plot sequence diagram of method on ax"""
    def add_rect_with_label(ax, x0, x1, label, ev_type):
        rect = Rectangle((x0, -1), x1-x0, 2, color=EVENT_COLORS[ev_type], alpha=0.4)
        ax.add_patch(rect)
        if label is not None:
            ax.annotate(label, ((x1-x0)/2, 0), color="black", ha="center", va="center")

    plot_mix = False
    idx = 0
    for i, row in df.iterrows():
        if row["type"] == "MixingEvent":
            plot_mix = True
            continue
        if plot_mix:
            # NOTE: What should labels be for mxixing events?
            add_rect_with_label(ax, x_data[idx], x_data[idx+1], "", "MixingEvent")
            idx += 1
        add_rect_with_label(ax, x_data[idx], x_data[idx+1], row["label"], row["type"])
        idx += 1
        plot_mix = False


def _plot_p_or_d(ax, x_data, y_data, name):
    """Helper method to plot p or d data on axis"""
    # Clean nan from y_data
    y_data = np.array(y_data[[np.any(~np.isnan(d)) for d in y_data]].tolist())
    # print(len(y_data), y_data)

    # Ensure no errors thrown when plotting data
    if len(x_data) - 1 != len(y_data) * 2:
        # Warn misshapen data?
        return

    for i in range(len(y_data[0])):
        ax.plot(
            x_data,
            [0] + [_ for num in y_data[:, i] for _ in (num, num)],
            color=COLORS[i],
            alpha=0.8,
            linewidth=2,
        )
    _range = max(  # Furthest point from origin
        abs(min([num for tup in y_data for num in tup])),
        abs(max([num for tup in y_data for num in tup])),
    )
    if name == "d" or _range > 5:
        # Even integers on y-axis
        _range = _range if _range % 2 == 0 else _range + 1  # Round up to next even int
        ax.set_ylim(-(_range + 2), _range + 2)
        ax.set_yticks(np.arange(-_range, _range + 1, step=2))
    else:
        ax.set_ylim(-(_range + 1), _range + 1)
        ax.set_yticks(np.arange(-_range, _range + 1, step=1))
    ax.grid(axis="y", color="black", alpha=0.2)
    ax.set_ylabel(name)


def _plot_data(ax, x_data, y_data, name):
    """Helper method to plot data and do formatting on ax"""
    # Clean nan from y_data
    y_data = np.array(y_data[[np.any(~np.isnan(d)) for d in y_data]].tolist())
    # print(len(y_data), y_data)

    # Ensure no errors thrown when plotting data
    if len(x_data) != len(y_data) * 2:
        # Warn misshapen data?
        return

    ax.plot(
        x_data,
        [_ for num in y_data for _ in (num, num)],
        color="b",
        alpha=0.6,
    )
    ax.set_ylabel(name.replace("_", " "))


def plot(df) -> plt.figure:
    """Create Plotly figure of symmetry pathways for method

    Args:
        DataFrame df: dataframe representation of method events to be plot

    Returns:
        figure fig: Matplotlib figure

    Example:
        TODO add example code
    """
    params = list(df.columns.drop(["type", "label", "duration", "fraction"]))
    # Move p and d to front of params
    if "d" in params:
        params.insert(0, params.pop(params.index("d")))
    if "p" in params:
        params.insert(0, params.pop(params.index("p")))

    if len(params) == 0:
        # Warn passed empty dataframe?
        return plt.figure()

    x_data = _make_x_data(df)
    offset_x_data = _offset_x_data(df, x_data)
    fig, axs = plt.subplots(
        nrows=len(params) + 1,
        ncols=1,
        figsize=(max(x_data) * 2.5, len(params) * 1.5 + 2),
        sharex=True,
        gridspec_kw={"hspace": 0.0},
    )

    axs[0].set_xlim(0, x_data[-1])
    axs[-1].get_xaxis().set_visible(False)

    # Plot sequence diagram
    _plot_sequence_diagram(axs[0], offset_x_data, df)

    # Iterate through axes and plot data
    for i, ax in enumerate(axs[1:], 0):
        if params[i] == "p" or params[i] == "d":
            _plot_p_or_d(ax, offset_x_data, df[params[i]], params[i])
        else:
            _plot_data(ax, x_data, df[params[i]], params[i])

    fig.tight_layout()
    return fig
