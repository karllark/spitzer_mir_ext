import argparse

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib
from astropy.table import Table

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

    avefilename = "test.fits"

    fontsize = 18

    font = {"size": fontsize}

    matplotlib.rc("font", **font)

    matplotlib.rc("lines", linewidth=1)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("xtick.minor", width=2)
    matplotlib.rc("ytick.major", width=2)
    matplotlib.rc("ytick.minor", width=2)

    figsize = (10, 8)
    fig, ax = pyplot.subplots(nrows=1, ncols=1, figsize=figsize)

    # Milky Way observed extinction from
    # Rieke & Lebofsky (1985)
    GL85_wave = np.array(
        [
            13.0,
            12.5,
            12.0,
            11.5,
            11.0,
            10.5,
            10.0,
            9.5,
            9.0,
            8.5,
            8.0,
            4.8,
            3.5,
            2.22,
            1.63,
            1.23,
        ]
    )
    GL85_axav = np.array(
        [
            0.027,
            0.030,
            0.037,
            0.047,
            0.060,
            0.074,
            0.083,
            0.087,
            0.074,
            0.043,
            0.020,
            0.023,
            0.058,
            0.112,
            0.175,
            0.282,
        ]
    )
    if args.elkejk:
        GL85_y, ji, ki = get_elkejk_from_alav(GL85_wave, GL85_axav)
    else:
        GL85_y = GL85_axav
    ax.plot(GL85_wave, GL85_y, "bo", label="GalCenter; Rieke & Lebofsky (1985)")

    # Chiar & Tielens (2005)
    # table in A(l)/A(K) units
    CT05 = Table.read(
        "data/pixie_dust_chiar_2005_modified.dat", format="ascii.commented_header"
    )
    CT05_wave = CT05["wave"]
    if args.elkejk:
        CT05_y, ji, ki = get_elkejk_from_alav(CT05_wave, 0.12 * CT05["local"])
    else:
        CT05_y = 0.12 * CT05["local"]
    ax.plot(CT05_wave, CT05_y, "co", label="GalCenter; Chiar & Tielens (2005)")

    # Zasowski et al. 2009 (supercedes Indebetouw?)
    # for l=90
    Z09_wave = np.array([1.22, 1.63, 2.22, 3.6, 4.5, 5.8, 8.0])
    Z09_ehlehk_l10 = np.array([-1.97, 0.0, 1.0, 1.79, 1.92, 2.08, 2.00])
    Z09_ehlehk_l90 = np.array([-2.18, 0.0, 1.0, 1.94, 2.38, 2.78, 2.49])
    if args.elkejk:
        Z09_y_l10 = get_elkejk_from_ehlehk(Z09_wave, Z09_ehlehk_l10)
        Z09_y_l90 = get_elkejk_from_ehlehk(Z09_wave, Z09_ehlehk_l90)
        ax.plot(
            Z09_wave,
            Z09_y_l10,
            "go",
            fillstyle="none",
            markersize=10,
            markeredgewidth=2.0,
            label="GalPlane (l10-15); Zasowski et al. (2009)",
        )
        ax.plot(
            Z09_wave,
            Z09_y_l90,
            "go",
            markersize=10,
            markeredgewidth=2.0,
            label="GalPlane (l90); Zasowski et al. (2009)",
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
    if args.elkejk:
        X16_y = -1 * X16_eklejk
        ax.plot(X16_wave, X16_y, "ro", label="GalPlane; Xue et al. (2016)")

    # Fritz et al. 2011
    F11 = Table.read("data/fritz11_galcenter.dat", format="ascii.commented_header")
    F11_wave = F11["wave"]
    if args.elkejk:
        F11_y, ji, ki = get_elkejk_from_alav(F11_wave, 0.12 * F11["ext"] / 2.62)
    else:
        F11_y = 0.12 * F11["ext"] / 2.62
    ax.plot(F11_wave, F11_y, "mo", label="GalCenter; Fritz et. al. (2011)")

    # get G19
    G19 = ExtData()
    G19.read(avefilename)
    G19_wave = G19.waves["BAND"]
    if args.elkejk:
        G19_y, ji, ki = get_elkejk_from_alav(G19_wave, G19.exts["BAND"])
    else:
        G19_y = G19.exts["BAND"]
    ax.plot(
        G19_wave,
        G19_y,
        "ko",
        markersize=12,
        markeredgewidth=2.0,
        label="Diffuse ISM; this work",
    )
    # plot the IRS data
    G19_IRS_y = G19.exts["IRS"]
    if args.elkejk:
        G19_IRS_y = G19_IRS_y / (G19.exts["BAND"][ji] - G19.exts["BAND"][ki])
        G19_akejk = G19.exts["BAND"][ki] / (G19.exts["BAND"][ji] - G19.exts["BAND"][ki])
        G19_IRS_y = G19_IRS_y - G19_akejk
    gindxs = np.where(G19.npts["IRS"] > 0)
    ax.plot(G19.waves["IRS"][gindxs], G19_IRS_y[gindxs], "k-")

    # finishing plot details
    if args.elkejk:
        ax.set_xlim(1.0, 40.0)
        ax.set_ylim(-1.0, 1.1)
        ax.set_ylabel(r"$E(\lambda - K)/E(J - K)$", fontsize=1.3 * fontsize)
    else:
        ax.set_ylim(0.0, 0.25)
        ax.set_ylabel(r"$A(\lambda)/A(V)$", fontsize=1.3 * fontsize)

    ax.set_yscale("linear")
    ax.set_xscale("log")

    ax.set_xlabel(r"$\lambda$ [$\mu m$]")

    ax.tick_params("both", length=10, width=2, which="major")
    ax.tick_params("both", length=5, width=1, which="minor")

    ax.legend(fontsize=fontsize)

    fig.tight_layout()

    save_str = "_litcomp"
    if args.png:
        fig.savefig(avefilename.replace(".fits", save_str + ".png"))
    elif args.pdf:
        fig.savefig(avefilename.replace(".fits", save_str + ".pdf"))
    else:
        pyplot.show()
