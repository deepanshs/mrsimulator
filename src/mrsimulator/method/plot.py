# -*- coding: utf-8 -*-
from itertools import groupby

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

__author__ = "Matthew D. Giammar"
__email__ = "giammar.7@buckeyemail.osu.edu"


DURATION_WIDTH = 0.6
MIXING_WIDTH = 0.05
# TODO: Ensure cannot run out of colors
COLORS = list(mcolors.TABLEAU_COLORS)
EVENT_COLORS = {  # TODO: add colors
    "ConstantDurationEvent": "orange",
    "SpectralEvent": "g",
    "MixingEvent": "b",
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
    # NOTE: Mixing event at end of sequence are ignored
    offset_x = np.array([0] + x_data)
    ev_groups = [(_type, sum(1 for _ in group)) for _type, group in groupby(df["type"])]
    # Remove MixingEvents from end of sequence
    if ev_groups[-1][0] == "MixingEvent":
        ev_groups.pop()

    for i, (_type, num) in enumerate(ev_groups):
        if i == 0 and _type == "MixingEvent":
            offset_x[1] += MIXING_WIDTH * num * 2
            continue
        if _type == "MixingEvent":
            offset_x[i] -= MIXING_WIDTH * num
            offset_x[i + 1] += MIXING_WIDTH * num

    return offset_x


def _add_rect_with_label(ax, x0, x1, label, ev_type):
    """Add a rectangle between x0 and x1 representing event"""
    rect_kwargs = {"color": EVENT_COLORS[ev_type], "alpha": 0.2}
    anno_kwargs = {
        "color": "black",
        "ha": "center",
        "va": "center",
        "rotation": 90 if ev_type == "MixingEvent" else 0,
    }
    rect = Rectangle((x0, -1), x1 - x0, 2, **rect_kwargs)
    ax.add_patch(rect)
    if label is not None:
        ax.annotate(label, ((x1 + x0) / 2, 0.5), **anno_kwargs)


def _plot_sequence_diagram(ax, x_data, df):
    """Helper method to plot sequence diagram of method on ax"""
    ev_groups = [(_type, sum(1 for _ in group)) for _type, group in groupby(df["type"])]
    if ev_groups[-1][0] == "MixingEvent":
        x_data = np.append(x_data, x_data[-1])
        x_data[-2] -= MIXING_WIDTH * ev_groups[-1][1] * 2

    df_idx = 0
    x_idx = 0
    for _type, num in ev_groups:
        if _type == "MixingEvent":
            # Create even spacing between MixingEvent rectangles within offset width
            tmp = np.linspace(x_data[x_idx], x_data[x_idx + 1], num=num + 1)
            for j in range(num):
                # Plot each MixingEvent
                _add_rect_with_label(
                    ax, tmp[j], tmp[j + 1], df["label"][df_idx + j], "MixingEvent"
                )
            x_idx += 1
        else:
            for j in range(num):
                # Increment x_idx if no mixing event between events
                if x_data[x_idx] == x_data[x_idx + 1]:
                    x_idx += 1
                _add_rect_with_label(
                    ax,
                    x_data[x_idx],
                    x_data[x_idx + 1],
                    df["label"][df_idx + j],
                    df["type"][df_idx + j],
                )
                x_idx += 1

        df_idx += num

    ax.axis("off")


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
    params = list(
        df.columns.drop(["type", "label", "duration", "fraction", "freq_contrib"])
    )
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

    return fig
