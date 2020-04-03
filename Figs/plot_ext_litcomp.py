import argparse

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib

# import astropy.units as u

# from dust_extinction.shapes import P92
from dust_extinction.averages import (
    RL85_MWAvg,
    I05_MWAvg,
    CT06_MWGC,
    CT06_MWLoc,
    F11_MWGC,
)

from measure_extinction.extdata import ExtData


def get_elkejk_from_alav(x, y):
    """
    Covert from A(l)/A(V) to A(l-K)/A(J-K)
    """
    kindx = np.argsort(np.absolute(x - 2.22))[0]
    jindx = np.argsort(np.absolute(x - 1.25))[0]
    alejk = y / (y[jindx] - y[kindx])
    elkejk = alejk - alejk[kindx]

    return (elkejk, jindx, kindx)


def get_elkejk_from_ehlehk(x, y):
    """
    Covert from E(H-l)/E(H-K) to A(l-K)/A(J-K)
    """
    kindx = np.argsort(np.absolute(x - 2.22))[0]
    # hindx = np.argsort(np.absolute(x - 1.63))[0]
    jindx = np.argsort(np.absolute(x - 1.25))[0]
    ehkejk = -1.0 / (y[jindx] - y[kindx])
    ehlejk = y * ehkejk
    elkejk = ehlejk[kindx] - ehlejk

    return elkejk


if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--elkejk", help="plot E(lambda-K)/A(J-K)", action="store_true")
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

    litmods = [RL85_MWAvg(), I05_MWAvg(), CT06_MWGC(), CT06_MWLoc(), F11_MWGC()]
    litdesc = [
        "GalCenter: Rieke & Lebofsky (1985)",
        "GalPlane: Indebetouw et al. (2005)",
        "GalCenter: Chiar & Tielens (2006)",
        "Local: Chiar & Tielens (2006)",
        "GalCenter: Fritz et al. (2011)",
    ]
    litfmt = ["bs", "gP", "c--", "c:", "m^"]
    for k, cmod in enumerate(litmods):
        lit_wave = 1.0 / cmod.obsdata_x
        lit_axav = cmod.obsdata_axav

        lit_y, ji, ki = get_elkejk_from_alav(lit_wave, lit_axav)
        ax[1].plot(
            lit_wave,
            lit_axav / lit_axav[ki],
            litfmt[k],
            alpha=0.25,
            label=litdesc[k],
            lw=3,
        )
        if k == 1:
            ax[0].plot(lit_wave, lit_y, litfmt[k], alpha=0.25, lw=3, label=litdesc[k])

    # Zasowski et al. 2009 (supercedes Indebetouw?)
    Z09_wave = np.array([1.22, 1.63, 2.22, 3.6, 4.5, 5.8, 8.0])
    Z09_ehlehk_l10 = np.array([-1.97, 0.0, 1.0, 1.79, 1.92, 2.08, 2.00])
    Z09_ehlehk_l90 = np.array([-2.18, 0.0, 1.0, 1.94, 2.38, 2.78, 2.49])

    Z09_y_l10 = get_elkejk_from_ehlehk(Z09_wave, Z09_ehlehk_l10)
    Z09_y_l90 = get_elkejk_from_ehlehk(Z09_wave, Z09_ehlehk_l90)
    ax[0].plot(
        Z09_wave,
        Z09_y_l10,
        "rs",
        fillstyle="none",
        markersize=8,
        markeredgewidth=2.0,
        alpha=0.5,
        label="GalPlane (l=10-15): Zasowski et al. (2009)",
    )
    ax[0].plot(
        Z09_wave,
        Z09_y_l90,
        "rs",
        markersize=8,
        markeredgewidth=1.0,
        alpha=0.5,
        label="GalPlane (l=90): Zasowski et al. (2009)",
    )

    # Xue et al. 2016
    X16_wave = np.array(
        [3.4, 4.6, 12.0, 22.0, 8.23, 3.6, 4.5, 5.8, 8.0, 24.0, 1.25, 2.22]
    )
    X16_eklejk = np.array(
        [
            0.238,
            0.312,
            0.269,
            0.370,
            0.273,
            0.260,
            0.313,
            0.355,
            0.334,
            0.428,
            -1.0,
            0.0,
        ]
    )
    X16_y = -1 * X16_eklejk
    ax[0].plot(X16_wave, X16_y, "yv", label="GalPlane: Xue et al. (2016)")

    # New measurements
    avefilename = "../data/all_ext_18feb20_ave.fits"

    # get G20_MWAvg
    G20 = ExtData()
    G20.read(avefilename)
    G20_wave = G20.waves["BAND"].value

    G20_ext = G20.exts["BAND"]
    G20_ext_uncs = G20.uncs["BAND"]

    gindxs_IRS = np.where(G20.npts["IRS"] > 0)
    G20_IRS_wave = G20.waves["IRS"][gindxs_IRS].value
    G20_IRS_ext = G20.exts["IRS"][gindxs_IRS]
    G20_IRS_uncs = G20.uncs["IRS"][gindxs_IRS]

    # photometry
    G20_y, ji, ki = get_elkejk_from_alav(G20_wave, G20_ext)

    ax[1].errorbar(
        G20_wave,
        G20_ext / G20_ext[ki],
        yerr=G20_ext_uncs / G20_ext[ki],
        fmt="ko",
        markersize=10,
        markeredgewidth=1.0,
        alpha=0.5,
        label="this work (JHK, IRAC, IRS15, MIPS24)",
    )

    ax[0].plot(
        G20_wave,
        G20_y,
        "ko",
        markersize=10,
        markeredgewidth=1.0,
        alpha=0.5,
        label="this work (JHK, IRAC, IRS15, MIPS24)",
    )

    # IRS
    ax[1].plot(
        G20_IRS_wave,
        G20_IRS_ext / G20_ext[ki],
        "k-",
        label="this work (IRS)",
        lw=2,
        alpha=0.65,
    )

    G20_IRS_y = G20_IRS_ext / (G20_ext[ji] - G20_ext[ki])
    G20_akejk = G20_ext[ki] / (G20_ext[ji] - G20_ext[ki])
    G20_IRS_y = G20_IRS_y - G20_akejk
    ax[0].plot(
        G20_IRS_wave, G20_IRS_y, "k-", lw=2, alpha=0.65, label="this work (IRS)",
    )

    for i in range(2):
        ax[i].set_yscale("linear")
        ax[i].set_xscale("log")
        ax[i].set_xlabel(r"$\lambda$ [$\mu m$]")
        ax[i].tick_params("both", length=10, width=2, which="major")
        ax[i].tick_params("both", length=5, width=1, which="minor")

        ax[i].legend(fontsize=0.75 * fontsize)

    # finishing plot details
    ax[0].set_xlim(1.0, 40.0)
    ax[0].set_ylim(-0.7, 1.1)
    ax[0].set_ylabel(r"$E(\lambda - K)/E(J - K)$", fontsize=1.3 * fontsize)

    ax[1].set_xlim(1.0, 40.0)
    ax[1].set_ylim(0.0, 3.0)
    ax[1].set_ylabel(r"$A(\lambda)/A(K)$", fontsize=1.3 * fontsize)
    ax[1].yaxis.tick_right()
    ax[1].yaxis.set_label_position("right")

    fig.tight_layout()

    save_str = "_litcomp"
    if args.png:
        fig.savefig(avefilename.replace(".fits", save_str + ".png"))
    elif args.pdf:
        fig.savefig(avefilename.replace(".fits", save_str + ".pdf"))
    else:
        pyplot.show()
