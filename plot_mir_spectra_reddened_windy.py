import argparse

import matplotlib.pyplot as plt
import matplotlib

import astropy.units as u

from plot_mir_spectra import plot_mir_set, ann_set

if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--eps", help="save figure as an eps file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    fontsize = 12

    font = {"size": fontsize}

    matplotlib.rc("font", **font)

    matplotlib.rc("lines", linewidth=1)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("xtick.minor", width=2)
    matplotlib.rc("ytick.major", width=2)
    matplotlib.rc("ytick.minor", width=2)

    # fig, ax = pyplot.subplots(nrows=1, ncols=1, figsize=(10, 13))
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))

    # kxrange = [3.0, 130.]
    kxrange = [3.0, 50.0]
    # ann_xvals = [41.0, 50.0]
    ann_xvals = [35.0, 42.0]
    spec_name = "IRS"
    norm_wave_range = [6.0, 10.0]
    # ann_wave_range = [15.0, 18.0]
    col_vals = ["b", "g", "r", "m", "c", "y"]
    # col_vals = ['k', 'k', 'k', 'k', 'k', 'k']

    # starnames = []
    # plot_mir_set(ax, starnames,
    #             ann_xvals=[6.0, 6.0], ann_wave_range=[5.0, 9.0])

    # ann_set(ax, fontsize, [35.0, 38.0], [3.9, 0.9], [45.0, 42.0],
    #        'Main Sequence', '(ordered by spectral type)')

    starnames = [
        "hd096042",
        "bd+63d1964",
        "hd149404",
        "hd169454",
        "hd229059",
        "hd166734",
        "vicyg5",
        "hd034921",
        "hd197702",
        "hd206773",
        "hd152408",
    ]
    plot_mir_set(
        ax,
        starnames,
        extra_off_val=0.0,
        ann_xvals=[9.0, 9.0],
        ann_wave_range=[6.0, 10.0] * u.micron,
    )

    ann_set(
        ax,
        fontsize,
        [35.0, 38.0],
        [10.5, 1.5],
        [45.0, 42.0],
        "Unusable for Extinction",
        "(ordered by wind srength)",
    )

    ax.set_yscale("linear")
    ax.set_ylim(0.5, 11.0)
    ax.set_xscale("log")
    ax.set_xlim(kxrange)
    ax.set_xlabel("$\lambda$ [$\mu m$]", fontsize=1.3 * fontsize)
    ax.set_ylabel("$\lambda^4 F(\lambda)/F(8 \mu m)$ + offset", fontsize=1.3 * fontsize)

    ax.tick_params("both", length=10, width=2, which="major")
    ax.tick_params("both", length=5, width=1, which="minor")

    fig.tight_layout()

    save_file = "spit_ir_mspec_reddened_windy"
    if args.png:
        fig.savefig(save_file + ".png")
    elif args.eps:
        fig.savefig(save_file + ".eps")
    elif args.pdf:
        fig.savefig(save_file + ".pdf")
    else:
        plt.show()
