# -*- coding: utf-8 -*-
import csdmpy as cp
import nmrglue as ng

__author__ = "Alexis McCarthy"
__email__ = "mccarthy.677@osu.edu"


def load_bruker_1d(filename):
    dic, data = ng.bruker.read(filename)
    udic = ng.bruker.guess_udic(dic, data)

    dims = list()
    for key, value in list(udic.items()):
        if type(key) == int:
            obj = cp.Dimension(
                type="linear",
                count=value["size"],
                increment=(str(1 / value["sw"]) + " s"),
                reciprocal={
                    "coordinates_offset": f'{value["car"]} Hz',
                    "origin_offset": f'{value["obs"]} MHz',
                },
                label=value["label"],
            )
            dims.append(obj)
    # dims = [cp.LinearDimension(
    #     count=value['size'],
    #     increment=(str(1 / value['sw']) + ' s'),
    #     coordinates_offset=(str(1 / value['car']) + ' s'),
    #     origin_offset=(str(1e-6 / value['obs']) + ' s'),
    #     label=value['label'])
    #     for key, value in list(udic.items()) if type(key) == int
    # ]
    csdm_data = cp.CSDM(
        dimensions=dims, dependent_variables=[cp.as_dependent_variable(data)]
    )
    return csdm_data, udic


# def plot_bruker_2d(data_2d):
#    x = data_2d.dimensions
#    y = data_2d.dependent_variables
#
#    x0 = x[0].coordinates
#    x1 = x[1].coordinates
#
#    y00 = y[0].components[0]
#
#    si = x[0].increment
#    extent = (
#        (x0[0] - 0.5 * si).value,
#        (x0[-1] + 0.5 * si).value,
#        x1[0].value,
#        x1[-1].value,
#    )
#
#    fig, axi = plt.subplots(
#        2, 2, gridspec_kw={"width_ratios": [4, 1], "height_ratios": [1, 4]}
#    )
#
#    ax = axi[1, 0]
#    im = NonUniformImage(ax, interpolation="nearest",
#                         extent=extent, cmap="bone_r")
#    im.set_data(x0, x1, y00.real / y00.real.max())
#
#    cbar = fig.colorbar(im)
#    cbar.ax.set_ylabel(y[0].axis_label[0])
#
#    ax.images.append(im)
#    for i in range(x1.size):
#        ax.plot(x0, np.ones(x0.size) * x1[i], "k--", linewidth=0.5)
#    ax.grid(axis="x", color="k", linestyle="--",
#            linewidth=0.5, which="both")
#
#    ax.set_xlim([extent[0], extent[1]])
#    ax.set_ylim([extent[2], extent[3]])
#    ax.set_xlabel(x[0].axis_label)
#    ax.set_ylabel(x[1].axis_label)
#    ax.set_title(y[0].name)
#
#    ax0 = axi[0, 0]
#    top = y00[-1].real
#    ax0.plot(x0, top, "k", linewidth=0.5)
#    ax0.set_xlim([extent[0], extent[1]])
#    ax0.set_ylim([top.min(), top.max()])
#    ax0.axis("off")
#
#    ax1 = axi[1, 1]
#    right = y00[:, 513].real
#    ax1.plot(right, x1, "k", linewidth=0.5)
#    ax1.set_ylim([extent[2], extent[3]])
#    ax1.set_xlim([right.min(), right.max()])
#    ax1.axis("off")
#
#    axi[0, 1].axis("off")
#
#    plt.tight_layout(pad=0.0, w_pad=0.0, h_pad=0.0)
#    plt.subplots_adjust(wspace=0.025, hspace=0.05)
#    plt.show()


# csdm_data, universal_dict = load_bruker_1d("ba2p2o7_zg")
# csdm_data, universal_dict = load_bruker_1d("bruker_1d")


# pp.pprint(universal_dict)
# print(csdm_data)
# plot_bruker_2d(data_2d)

# ax = plt.subplot(projection="csdm")
# ax.plot(csdm_data)
# plt.tight_layout()
# plt.show()


# x = range(len(csdm_data))
# plt.plot(csdm_data)
# plt.show()
