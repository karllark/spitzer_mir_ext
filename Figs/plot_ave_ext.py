import argparse

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib
from matplotlib.ticker import ScalarFormatter

import astropy.units as u

from dust_extinction.shapes import P92
from dust_extinction.parameter_averages import F19

from measure_extinction.extdata import ExtData


if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    fontsize = 14

    font = {"size": fontsize}

    matplotlib.rc("font", **font)

    matplotlib.rc("lines", linewidth=1)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("xtick.minor", width=2)
    matplotlib.rc("ytick.major", width=2)
    matplotlib.rc("ytick.minor", width=2)

    figsize = (14, 6)
    fig, ax = pyplot.subplots(nrows=1, ncols=2, figsize=figsize)

    # F19 R(V) average model
    mod_x = np.logspace(np.log10(0.115), np.log10(2.5), num=1000) * u.micron
    F19_Rv = F19(Rv=3.1)
    for cax in ax:
        cax.plot(mod_x, F19_Rv(mod_x), "k--", lw=2, alpha=0.65, label="F19 R(V) = 3.1")

    # New measurements
    avefilenames = [
        "data/all_ext_18feb20_diffuse_ave.fits",
        # "fits/hd283809_hd064802_ext_P92_FM90.fits",
        # "fits/hd029647_hd195986_ext_P92_FM90.fits",
    ]
    pcol = ["b", "g", "c"]
    psym = ["o", "s", "^"]
    pline = ["-", ":", "-."]
    plabel = ["diffuse", "HD283809", "HD029647"]

    for i, avefilename in enumerate(avefilenames):
        G20 = ExtData()
        G20.read(avefilename)
        G20_wave = G20.waves["BAND"].value

        G20_ext = G20.exts["BAND"]
        G20_ext_uncs = G20.uncs["BAND"]

        gindxs_IRS = np.where(G20.npts["IRS"] > 0)
        G20_IRS_wave = G20.waves["IRS"][gindxs_IRS].value
        G20_IRS_ext = G20.exts["IRS"][gindxs_IRS]
        G20_IRS_uncs = G20.uncs["IRS"][gindxs_IRS]

        gindxs_IUE = np.where(G20.npts["IUE"] > 0)
        G20_IUE_wave = G20.waves["IUE"][gindxs_IUE].value
        G20_IUE_ext = G20.exts["IUE"][gindxs_IUE]
        G20_IUE_uncs = G20.uncs["IUE"][gindxs_IUE]

        ax[1].errorbar(
            G20_wave,
            G20_ext,
            yerr=G20_ext_uncs,
            fmt=pcol[i] + psym[i],
            markersize=10,
            markeredgewidth=1.0,
            alpha=0.5,
            label=plabel[i],  # + " (JHK, IRAC, IRS15, MIPS24)",
        )

        # IRS
        ax[1].plot(
            G20_IRS_wave,
            G20_IRS_ext,
            pcol[i] + pline[i],
            # label=plabel[i] + " (IRS)",
            lw=2,
            alpha=0.65,
        )

        ax[0].errorbar(
            G20_wave,
            G20_ext,
            yerr=G20_ext_uncs,
            fmt=pcol[i] + psym[i],
            markersize=10,
            markeredgewidth=1.0,
            alpha=0.5,
            label=plabel[i],  # + " (UBV)",
        )

        # IUE
        ax[0].plot(
            G20_IUE_wave,
            G20_IUE_ext,
            pcol[i] + pline[i],
            # label=plabel[i] + " (IUE/STIS)",
            lw=2,
            alpha=0.85,
        )

    for i in range(2):
        ax[i].set_yscale("linear")
        ax[i].set_xlabel(r"$\lambda$ [$\mu m$]")
        ax[i].tick_params("both", length=10, width=2, which="major")
        ax[i].tick_params("both", length=5, width=1, which="minor")
        ax[i].set_ylabel(r"$A(\lambda)/A(V)$", fontsize=1.3 * fontsize)

        ax[i].legend(fontsize=0.75 * fontsize)

    # finishing plot details
    ax[0].set_xlim(0.1, 0.6)
    ax[0].set_xscale("log")
    ax[0].set_ylim(0.0, 4.0)
    ax[0].xaxis.set_major_formatter(ScalarFormatter())
    ax[0].xaxis.set_minor_formatter(ScalarFormatter())

    ax[1].set_xscale("log")
    ax[1].set_xlim(1.0, 40.0)
    ax[1].set_ylim(0.0, 0.4)
    ax[1].yaxis.tick_right()
    ax[1].yaxis.set_label_position("right")
    ax[1].xaxis.set_major_formatter(ScalarFormatter())

    fig.tight_layout()

    save_fname = "diffuse_aveext"
    if args.png:
        fig.savefig(save_fname + ".png")
    elif args.pdf:
        fig.savefig(save_fname + ".pdf")
    else:
        pyplot.show()
