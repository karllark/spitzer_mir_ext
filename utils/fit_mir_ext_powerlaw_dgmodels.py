import matplotlib.pyplot as plt
import matplotlib

import numpy as np
import argparse
import warnings

import astropy.units as u

from astropy.modeling.fitting import LevMarLSQFitter

from G21 import G21, G21_drude_asym

from dust_extinction.grain_models import D03, ZDA04, J13


def clean_pnames(pnames):
    """
    function to clean of the _? part of the names due to making a CompoundModel
    """
    if pnames[0][-1] in ["0", "1"]:
        clean_pnames = [cpname[:-2] for cpname in pnames]
        return clean_pnames
    else:
        return pnames


if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dgmodel",
        help="dust grain model to use",
        choices=["D03", "ZDA04", "J13"],
        default="D03",
    )
    parser.add_argument(
        "--symfit", help="include symmetric drude fit", action="store_true"
    )
    parser.add_argument("--notitle", help="no title on plot", action="store_true")
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    parser.add_argument("--path", help="path for the extinction curves")
    args = parser.parse_args()

    wave = np.logspace(np.log10(1.0), np.log10(39.0), num=100) * u.micron
    if args.dgmodel == "D03":
        dmod = D03()
    elif args.dgmodel == "ZDA04":
        dmod = ZDA04()
    elif args.dgmodel == "J13":
        dmod = J13()
    else:
        print(f"{args.dgmodel} model not recognized")
        exit()

    y = dmod(wave)
    ofile = "fits/D03.fits"

    # remove units as fitting routines often cannot take numbers with units
    x = wave.to(1.0 / u.micron, equivalencies=u.spectral()).value

    g21_init = G21()
    g21_asym_init = G21_drude_asym()

    # fit the extinction only using data between 1 and 40 micron
    gvals = (1.0 < 1.0 / x) & (1.0 / x < 40.0)
    fit = LevMarLSQFitter()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        g21_fit = fit(g21_init, x[gvals], y[gvals],)

        g21_asym_fit = fit(g21_asym_init, x[gvals], y[gvals],)

    # save the extinction curve and fit
    best_params = (clean_pnames(g21_asym_fit.param_names), g21_asym_fit.parameters)

    print(best_params)

    # setup the plot
    fontsize = 18
    # fontsize = 10
    font = {"size": fontsize}
    matplotlib.rc("font", **font)
    matplotlib.rc("lines", linewidth=1.5)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("xtick.minor", width=2)
    matplotlib.rc("ytick.major", width=2)
    matplotlib.rc("ytick.minor", width=2)

    fig, ax = plt.subplots(
        nrows=2, figsize=(12, 8), sharex=True, gridspec_kw={"height_ratios": [5, 1]}
    )

    ax[0].plot(wave, y, "k-")
    # obsext.plot(ax[0], color="k")
    g21_fit_y = g21_fit(wave[gvals])
    g21_asym_fit_y = g21_asym_fit(wave[gvals])

    ax[0].set_ylabel(r"$A(\lambda)/A(V)$", fontsize=1.3 * fontsize)

    ax[1].set_xlabel(r"$\lambda$ [$\mu m$]", fontsize=1.3 * fontsize)

    if args.symfit:
        ax[0].plot(wave[gvals], g21_fit_y, "g-", label="Sym Best Fit", alpha=0.7)
    ax[0].plot(wave[gvals], g21_asym_fit_y, "b-", label="Best Fit")

    mmy = np.array([min(g21_fit_y), max(g21_fit_y)])
    mmd = 0.1 * (mmy[1] - mmy[0])
    ax[0].set_ylim(mmy + np.array([-1.0, 1.0]) * mmd)
    ax[0].set_xlim(1.0, 40.0)
    ax[0].set_xscale("log")
    if not args.notitle:
        ax[0].set_title(ofile)

    g21_comps = g21_asym_fit.copy()
    g21_comps.sil1_amp = 0.0
    ax[0].plot(wave[gvals], g21_comps(wave[gvals]), "k--", alpha=0.5)

    g21_comps = g21_asym_fit.copy()
    g21_comps.sil2_amp = 0.0
    ax[0].plot(wave[gvals], g21_comps(wave[gvals]), "k--", alpha=0.5)

    g21_comps = g21_asym_fit.copy()
    g21_comps.sil1_amp = 0.0
    g21_comps.sil2_amp = 0.0
    ax[0].plot(wave[gvals], g21_comps(wave[gvals]), "k--", alpha=0.5)

    ax[0].legend(loc="best")

    # residuals
    ax[1].plot(wave[gvals], np.zeros((len(wave[gvals]))), "k--")
    if args.symfit:
        ax[1].plot(wave[gvals], y[gvals] - g21_fit_y, "g-", alpha=0.7)

    ax[1].plot(
        wave, y - g21_asym_fit(wave), "b-",
    )
    ax[1].set_ylim(np.array([-1.0, 1.0]) * mmd)

    plt.tight_layout()

    # plot or save to a file
    outname = ofile.replace(".fits", "")
    if args.png:
        fig.savefig(outname + ".png")
    elif args.pdf:
        fig.savefig(outname + ".pdf")
    else:
        plt.show()
