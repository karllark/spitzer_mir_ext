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

    red_starnames = [
        "hd281159",
        "hd283809",
        "hd029309",
        "hd147889",
        "hd204827",
        "vicyg1",
        "hd029647",
        "hd147701",
        "hd014956",
        "vicyg2",
        "hd112272",
        "bd+63d1964",
        "hd192660",
        "hd229238",
        "vicyg8a",
    ]

    windy_starnames = [
        "hd096042",
        "hd147933",
        "hd197702",
        "hd149404",
        "hd169454",
        "hd229059",
        "hd166734",
        "hd206773",
        "vicyg5",
        "hd034921",
        "hd152408",
    ]

    path = "/home/kgordon/Python_git/extstar_data/"
    subpath = "DAT_files/"

    xbands = ("J", "K")
    ybands = ("K", "MIPS24")

    # standard stars
    xvals_comp, yvals_comp = get_colors(comp_starnames)
    xvals_red, yvals_red = get_colors(red_starnames)
    xvals_windy, yvals_windy = get_colors(windy_starnames)

    fontsize = 18

    font = {"size": fontsize}

    plt.rc("font", **font)

    plt.rc("lines", linewidth=1)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("xtick.minor", width=2)
    plt.rc("ytick.major", width=2)
    plt.rc("ytick.minor", width=2)

    fig, ax = plt.subplots(figsize=(10, 8))

    ax.errorbar(
        xvals_comp[:, 0],
        yvals_comp[:, 0],
        xerr=xvals_comp[:, 1],
        yerr=yvals_comp[:, 1],
        fmt="ko",
        label="unreddened",
        alpha=0.5,
    )
    ax.errorbar(
        xvals_red[:, 0],
        yvals_red[:, 0],
        xerr=xvals_red[:, 1],
        yerr=yvals_red[:, 1],
        fmt="go",
        label="reddened",
    )
    ax.errorbar(
        xvals_windy[:, 0],
        yvals_windy[:, 0],
        xerr=xvals_windy[:, 1],
        yerr=yvals_windy[:, 1],
        fmt="ro",
        label="windy+bad",
        alpha=0.5,
    )

    ax.set_xlabel(f"{xbands[0]} - {xbands[1]}")
    ax.set_ylabel(f"{ybands[0]} - {ybands[1]}")

    ax.legend()

    fig.tight_layout()

    basefile = "irext_nir_mir_phot"
    if args.png:
        fig.savefig(basefile + ".png")
    elif args.pdf:
        fig.savefig(basefile + ".pdf")
    else:
        plt.show()
