#!/usr/bin/env python
#
# Program to plot a list of extinction curves
#
from __future__ import absolute_import, division, print_function, unicode_literals
import argparse

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib

import astropy.units as u

from measure_extinction.extdata import ExtData


def get_elkejk_from_elv(x_band, elv_band, x, elv):
    """
    Covert from E(l-V) to E(l-K)/E(J-K)
    """
    kindx = np.argsort(np.absolute(x_band - 2.22))[0]
    jindx = np.argsort(np.absolute(x_band - 1.25))[0]
    elk = elv - elv_band[kindx]
    elkejk = elk / (elv_band[jindx] - elv_band[kindx])

    return elkejk


if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("filelist", help="file with list of curves to plot")
    parser.add_argument(
        "--rebin_fac", type=int, default=None, help="rebin factor for spectra"
    )
    parser.add_argument(
        "--prevobs", help="plot previous observations", action="store_true"
    )
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    filename = args.filelist

    f = open(filename, "r")
    file_lines = list(f)
    extnames = []
    extdatas = []
    avs = []
    for line in file_lines:
        if (line.find("#") != 0) & (len(line) > 0):
            name = line.rstrip()
            extnames.append(name)
            text = ExtData(filename="fits/%s" % name)
            extdatas.append(text)
            avs.append(text.columns["AV"][0])

    fontsize = 18

    font = {"size": fontsize}

    matplotlib.rc("font", **font)

    matplotlib.rc("lines", linewidth=1)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("xtick.minor", width=2)
    matplotlib.rc("ytick.major", width=2)
    matplotlib.rc("ytick.minor", width=2)

    figsize = (10, 12)
    fig, ax = pyplot.subplots(nrows=1, ncols=1, figsize=figsize)

    sindxs = np.argsort(avs)

    col_vals = ["b", "g", "r", "m", "c", "y"]
    lin_vals = ["--", ":", "-."]

    mod_x = 1.0 / np.arange(1.0, 40.0, 0.1)
    for i in range(len(extnames)):
        k = sindxs[i]

        color = col_vals[i % 6]

        curtype = "BAND"
        gindxs, = np.where(extdatas[k].npts[curtype] > 0)
        x_band = extdatas[k].waves[curtype][gindxs].to(u.micron).value
        elv_band = extdatas[k].exts[curtype][gindxs]

        for curtype in extdatas[k].waves.keys():
            gindxs, = np.where(extdatas[k].npts[curtype] > 0)
            x = extdatas[k].waves[curtype][gindxs].to(u.micron).value
            y = extdatas[k].exts[curtype][gindxs]
            yu = extdatas[k].uncs[curtype][gindxs]
            y_new = get_elkejk_from_elv(x_band, elv_band, x, y)

            if len(gindxs) < 20:
                # plot small number of points (usually BANDS data) as
                # points
                ax.plot(x, y_new, "o", color=color, mfc="white")
            else:
                ax.plot(x, y_new, "-", color=color)

    ax.set_yscale("linear")
    ax.set_xscale("log")

    # ax.set_xlim(1.0, 100.0)
    # ax.set_ylim(-6, -1.0)
    ax.set_ylabel(r"$E(\lambda - K)/E(J-K)$", fontsize=1.3 * fontsize)
    ax.set_xlabel(r"$\lambda$ [$\mu m$]")

    ax.tick_params("both", length=10, width=2, which="major")
    ax.tick_params("both", length=5, width=1, which="minor")

    fig.tight_layout()

    save_str = "_mext_elkejk"
    if args.png:
        fig.savefig(args.filelist.replace(".dat", save_str + ".png"))
    elif args.pdf:
        fig.savefig(args.filelist.replace(".dat", save_str + ".pdf"))
    else:
        pyplot.show()
