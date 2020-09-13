import argparse
import warnings
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

import astropy.units as u
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.modeling.models import Polynomial1D, Drude1D

# from dust_extinction.conversions import AxAvToExv
from measure_extinction.extdata import ExtData

# from models_mcmc_extension import EmceeFitter


if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("extfile", help="file with extinction curve")
    parser.add_argument(
        "--burnfrac", type=float, default=0.1, help="fraction of MCMC chain to burn"
    )
    parser.add_argument(
        "--nsteps", type=int, default=100, help="# of steps in MCMC chain"
    )
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    # get a saved extnction curve
    file = args.extfile
    # file = '/home/kgordon/Python_git/spitzer_mir_ext/fits/hd147889_hd064802_ext.fits'
    ofile = file.replace(".fits", "_POLYDRUDE.fits")
    extdata = ExtData(filename=file)

    # get an observed extinction curve to fit
    (wave, y, y_unc) = extdata.get_fitdata(
        ["BAND", "IRS"], remove_uvwind_region=True, remove_lya_region=True
    )
    # remove units as fitting routines often cannot take numbers with units
    x = wave.to(1.0 / u.micron, equivalencies=u.spectral()).value

    # determine the initial guess at the A(V) values
    #  just use the average at wavelengths > 5
    #  limit as lambda -> inf, E(lamda-V) -> -A(V)
    (indxs,) = np.where(1.0 / x > 5.0)

    # initialize the model
    #    a few tweaks to the starting parameters helps find the solution
    ponly = (
        Polynomial1D(degree=5, c0=0.0)  # , fixed={"c0": True})
        + Drude1D(
            amplitude=0.1,
            x_0=0.1,
            fwhm=0.01,
            bounds={
                "amplitude": [0.0, None],
                "x_0": [1.0 / 10.5, 1.0 / 9.0],
                "fwhm": [0.0001, 0.5],
            },
        )
        + Drude1D(
            amplitude=0.01,
            x_0=1.0 / 20.0,
            fwhm=0.05,
            bounds={
                "amplitude": [0.0, None],
                "x_0": [1.0 / 22.0, 1.0 / 17.0],
                "fwhm": [0.0001, 0.5],
            },
        )
        #    + Gaussian1D(amplitude=1.0, mean=4.6, stddev=1.0,
        #        bounds={"amplitude": [0.0, 10.0], "mean": [4.5, 4.7], "stddev": [0.5, 1.5]})
    )

    ponly[0].c0 = 0.001
    ponly[0].c0.bounds = (0.0, None)
    ponly[0].c1 = 0.05
    ponly[0].c2 = -0.1
    ponly[0].c3 = 0.5
    ponly[0].c4 = -0.2
    ponly[0].c5 = 0.03

    # pick the fitter
    fit = LevMarLSQFitter()

    # fit the data to the P92 model using the fitter
    # p92_fit = fit(p92_init, x, y, weights=1.0 / y_unc, maxiter=1000)
    (gindxs,) = np.where(1.0 / x < 40.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        p92_fit = fit(ponly, x[gindxs], y[gindxs], weights=1.0 / y_unc[gindxs], maxiter=1000, epsilon=0.001)

    for k, cur_pname in enumerate(p92_fit.param_names):
        print("{:12} {:6.4e}".format(cur_pname, p92_fit.parameters[k]))

    # plotting setup for easier to read plots
    fontsize = 18
    font = {"size": fontsize}
    matplotlib.rc("font", **font)
    matplotlib.rc("lines", linewidth=1)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("xtick.minor", width=2)
    matplotlib.rc("ytick.major", width=2)
    matplotlib.rc("ytick.minor", width=2)

    # setup the plot
    fig, ax = plt.subplots(figsize=(12, 8))

    # subplot
    ax2 = plt.axes([0.60, 0.35, 0.35, 0.35])

    # plot the bands and all spectra for this star
    extdata.plot(ax, color="k", alpha=0.5)
    extdata.plot(ax2, color="k", alpha=0.5)

    # ax.plot(1.0 / x, p92_init(x), "r--", label="P92 Init")
    ax.plot(1.0 / x, ponly(x), "b-", label="initial guess")
    ax.plot(1.0 / x, p92_fit(x), "r-", label="Best Fit")
    ax2.plot(1.0 / x, p92_fit(x), "r-")

    # finish configuring the plot
    ax.set_yscale("linear")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\lambda$ [$\mu m$]", fontsize=1.3 * fontsize)
    ax.set_ylabel(extdata._get_ext_ytitle(extdata.type), fontsize=1.3 * fontsize)
    ax.tick_params("both", length=10, width=2, which="major")
    ax.tick_params("both", length=5, width=1, which="minor")
    ax.legend()

    # finish configuring the subplot
    sp_xlim = [1.0, 35.0]
    ax2.set_xlim(sp_xlim)
    # ax2.set_ylim(-best_fit_Av-0.1, -best_fit_Av+0.5)
    (indxs,) = np.where((x > 1.0 / sp_xlim[1]) & (x < 1.0 / sp_xlim[0]))
    # ax2.set_ylim(
    #    min([min(p92_fit(x)[indxs]), -best_fit_Av]) - 0.1, max(p92_fit(x)[indxs]) + 0.1
    # )
    ax2.set_yscale("log")
    ax2.set_ylim(0.01, 0.4)

    # use the whitespace better
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        fig.tight_layout()

    # plot or save to a file
    outname = ofile.replace(".fits", "")
    if args.png:
        fig.savefig(outname + ".png")
    elif args.pdf:
        fig.savefig(outname + ".pdf")
    else:
        plt.show()
