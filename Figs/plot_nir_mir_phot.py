import argparse
import numpy as np

import matplotlib.pyplot as plt
from measure_extinction.stardata import StarData


def get_colors(starnames):

    # standard stars
    n_stars = len(starnames)
    xvals = np.zeros((n_stars, 2))
    yvals = np.zeros((n_stars, 2))
    for i in range(n_stars):
        stardata = StarData(subpath + starnames[i] + ".dat", path=path, use_corfac=True)

        b1vals = stardata.data["BAND"].get_band_mag(xbands[0])
        b2vals = stardata.data["BAND"].get_band_mag(xbands[1])
        if (b1vals is not None) and (b2vals is not None):
            xvals[i, 0] = b1vals[0] - b2vals[0]
            xvals[i, 1] = np.sqrt((b1vals[1] ** 2) + (b2vals[1] ** 2))

        b1vals = stardata.data["BAND"].get_band_mag(ybands[0])
        b2vals = stardata.data["BAND"].get_band_mag(ybands[1])

        if (b1vals is not None) and (b2vals is not None):
            yvals[i, 0] = b1vals[0] - b2vals[0]
            yvals[i, 1] = np.sqrt((b1vals[1] ** 2) + (b2vals[1] ** 2))

    return (xvals, yvals)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--mirband", default="MIPS24", help="MIR band for color")
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    comp_starnames = [
        "hd064802",
        "hd074273",
        "hd031726",
        "hd036512",
        "hd214680",
        "hd047839",
        "hd195986",
        "hd034816",
        "hd165024",
        "hd064760",
        "hd204172",
        "hd188209",
        "hd051283",
    ]

    red_starnames_super = [
        "hd029647",
        "hd147701",
        "hd014956",
        "cygob2-2",
        "hd112272",
        "bd+63d1964",
        "hd192660",
        "hd229238",
        "cygob2-8a",
    ]

    red_starnames_main = [
        "hd147933",
        "hd281159",
        "hd283809",
        "hd029309",
        "hd147889",
        "hd204827",
        "cygob2-1",
    ]

    windy_starnames = [
        "hd197702",
        "hd149404",
        "hd169454",
        "hd229059",
        "hd166734",
        "hd206773",
        "cygob2-5",
        "hd034921",
        "hd152408",
    ]

    bad_starnames = ["hd096042"]

    path = "/home/kgordon/Python_git/extstar_data/"
    subpath = "DAT_files/"

    xbands = ("J", "K")
    ybands = ("K", args.mirband)
    lab_xbands = ("$J$", "$K_s$")
    lab_ybands = ("$K_s$", f"${args.mirband}$")

    # standard stars
    xvals_comp, yvals_comp = get_colors(comp_starnames)
    xvals_red_super, yvals_red_super = get_colors(red_starnames_super)
    xvals_red_main, yvals_red_main = get_colors(red_starnames_main)
    xvals_windy, yvals_windy = get_colors(windy_starnames)
    xvals_bad, yvals_bad = get_colors(bad_starnames)

    fontsize = 18

    font = {"size": fontsize}

    plt.rc("font", **font)

    plt.rc("lines", linewidth=1)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("xtick.minor", width=2)
    plt.rc("ytick.major", width=2)
    plt.rc("ytick.minor", width=2)

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.errorbar(
        xvals_comp[:, 0],
        yvals_comp[:, 0],
        xerr=xvals_comp[:, 1],
        yerr=yvals_comp[:, 1],
        fmt="ko",
        label="comparison",
        alpha=0.5,
    )
    ax.errorbar(
        xvals_red_super[:, 0],
        yvals_red_super[:, 0],
        xerr=xvals_red_super[:, 1],
        yerr=yvals_red_super[:, 1],
        fmt="go",
        label="reddened",
    )
    ax.errorbar(
        xvals_red_main[:, 0],
        yvals_red_main[:, 0],
        xerr=xvals_red_main[:, 1],
        yerr=yvals_red_main[:, 1],
        fmt="go",
        # label="reddened V-IV",
        # fillstyle="none",
    )
    ax.errorbar(
        xvals_windy[:, 0],
        yvals_windy[:, 0],
        xerr=xvals_windy[:, 1],
        yerr=yvals_windy[:, 1],
        fmt="ro",
        label="windy",
        alpha=0.5,
    )
    # ax.errorbar(
    #    xvals_bad[:, 0],
    #    yvals_bad[:, 0],
    #    xerr=xvals_bad[:, 1],
    #    yerr=yvals_bad[:, 1],
    #    fmt="co",
    #    label="bad",
    #    alpha=0.5,
    # )

    # line at 0.6
    if args.mirband == "MIPS24":
        ax.plot([-0.3, 0.9], [0.6, 0.6], "k--", alpha=0.25, lw=2)
    if args.mirband == "IRAC4":
        ax.plot([-0.3, 0.9], [0.4, 0.4], "k--", alpha=0.25, lw=2)

    ax.set_xlabel(f"{lab_xbands[0]} - {lab_xbands[1]}")
    ax.set_ylabel(f"{lab_ybands[0]} - {lab_ybands[1]}")

    ax.legend()

    fig.tight_layout()

    basefile = "irext_nir_mir_phot"
    if args.png:
        fig.savefig(basefile + ".png")
    elif args.pdf:
        fig.savefig(basefile + ".pdf")
    else:
        plt.show()
