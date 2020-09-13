import matplotlib.pyplot as plt
import matplotlib

import numpy as np
import argparse
import warnings

import astropy.units as u

# from astropy.modeling.models import PowerLaw1D, Drude1D
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.modeling.fitting import _fitter_to_model_params

from measure_extinction.extdata import ExtData
from dust_extinction.conversions import AxAvToExv

from models_mcmc_extension import EmceeFitter

from G21 import G21, G21_drude_asym


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
    parser.add_argument("file", help="file with the extinction curve to fit")
    parser.add_argument(
        "--burnfrac", type=float, default=0.1, help="fraction of MCMC chain to burn"
    )
    parser.add_argument(
        "--nsteps", type=int, default=100, help="# of steps in MCMC chain"
    )
    parser.add_argument(
        "--symfit", help="include symmetric drude fit", action="store_true"
    )
    parser.add_argument(
        "--notitle", help="no title on plot", action="store_true"
    )
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    parser.add_argument("--path", help="path for the extinction curves")
    args = parser.parse_args()

    if args.path:
        locpath = args.path + "/"
    else:
        locpath = ""

    file = args.file
    ofile = file.replace(".fits", "_POWLAW2DRUDE.fits")

    # read in the observed E(l-V) or A(l)/A(V) extinction curve
    obsext = ExtData(filename=locpath + file)

    # get an observed extinction curve to fit
    (wave, y, y_unc) = obsext.get_fitdata(["BAND", "IRS"])

    # remove units as fitting routines often cannot take numbers with units
    x = wave.to(1.0 / u.micron, equivalencies=u.spectral()).value

    if obsext.type == "elx":
        # determine the initial guess at the A(V) values
        #  just use the average at wavelengths > 5
        #  limit as lambda -> inf, E(lamda-V) -> -A(V)
        (indxs,) = np.where(1.0 / x > 5.0)
        av_guess = -1.0 * np.average(y[indxs])
        if not np.isfinite(av_guess):
            av_guess = 1.0

        g21_init = G21() | AxAvToExv(Av=av_guess)
        g21_asym_init = G21_drude_asym() | AxAvToExv(Av=av_guess)
        # g21_init = g21_x() | AxAvToExv(Av=av_guess)

        g21_asym_init[0].sil2_fwhm.fixed = True

        if max(1.0 / x) <= 20.0:
            print("fix amp")
            g21_asym_init[0].sil2_amp.fixed = True
            g21_asym_init[0].sil2_center.fixed = True
            g21_asym_init[0].sil2_fwhm.fixed = True
            g21_asym_init[0].sil2_asym.fixed = True
        elif max(1.0 / x) <= 26.0:
            # g21_asym_init[0].sil2_amp.fixed = True
            g21_asym_init[0].sil2_center.fixed = True
            g21_asym_init[0].sil2_fwhm.fixed = True
            g21_asym_init[0].sil2_asym.fixed = True
    elif obsext.type == "alav":
        g21_init = G21()
        g21_asym_init = G21_drude_asym()

        # g21_asym_init.sil2_fwhm.fixed = True

    # g21_asym_init[0].sil1_amp.fixed = True

    # fit the extinction only using data between 1 and 40 micron
    gvals = (1.0 < 1.0 / x) & (1.0 / x < 40.0)
    fit = LevMarLSQFitter()

    nsteps = args.nsteps
    emcee_samples_file = ofile.replace(".fits", ".h5")
    fit2 = EmceeFitter(
        nsteps=nsteps, burnfrac=args.burnfrac, save_samples=emcee_samples_file
    )

    weights = 1.0 / y_unc[gvals]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        g21_fit = fit(
            g21_init,
            x[gvals],
            y[gvals],
            weights=weights,
            # maxiter=10000,
            # epsilon=0.001,
        )

        g21_asym_fit = fit(
            g21_asym_init,
            x[gvals],
            y[gvals],
            weights=weights,
            # maxiter=10000,
            # epsilon=0.1,
        )
        # print(g21_fit.parameters)
        # print(g21_asym_fit.param_names)
        # print(g21_asym_fit.parameters)

        g21_asym_fit2 = fit2(g21_asym_fit, x[gvals], y[gvals], weights=weights)

    # make the standard mcmc plots
    fit2.plot_emcee_results(g21_asym_fit2, filebase=ofile.replace(".fits", ""))

    # save the extinction curve and fit
    best_params = (clean_pnames(g21_asym_fit.param_names), g21_asym_fit.parameters)
    per_param_vals = zip(
        g21_asym_fit2.parameters, g21_asym_fit2.uncs_plus, g21_asym_fit2.uncs_minus
    )
    per_params = (clean_pnames(g21_asym_fit2.param_names), list(per_param_vals))

    with warnings.catch_warnings():
        # warnings.simplefilter("ignore", category=AstropyWarning)
        g21_params = {"type": "G21", "best": best_params, "per": per_params}
        obsext.save(ofile, save_params=g21_params)

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

    obsext.plot(ax[0], color="k")
    g21_fit_y = g21_fit(wave[gvals])
    g21_asym_fit_y = g21_asym_fit(wave[gvals])

    if obsext.type == "elx":
        if args.symfit:
            ax[0].plot(
                wave[gvals],
                g21_fit.Av_1.value * np.full((len(wave[gvals])), -1.0),
                "g:",
                label="-A(V)",
            )
        ax[0].plot(
            wave[gvals],
            g21_asym_fit.Av_1.value * np.full((len(wave[gvals])), -1.0),
            "b:",
            label="-A(V)",
        )
        ax[0].set_ylabel(r"$E(\lambda - V)$", fontsize=1.3 * fontsize)
    else:
        ax[0].set_ylabel(r"$A(\lambda)/A(V)$", fontsize=1.3 * fontsize)

    ax[1].set_xlabel(r"$\lambda$ [$\mu m$]", fontsize=1.3 * fontsize)

    # plot samples from the mcmc chaing
    flat_samples = fit2.fit_info["sampler"].get_chain(
        discard=int(0.1 * nsteps), flat=True
    )
    inds = np.random.randint(len(flat_samples), size=100)
    model_copy = g21_asym_fit2.copy()
    for ind in inds:
        sample = flat_samples[ind]
        _fitter_to_model_params(model_copy, sample)
        ax[0].plot(wave[gvals], model_copy(wave[gvals]), "C1", alpha=0.05)
    # for the figure legend
    ax[0].plot(wave[gvals], g21_asym_fit2(wave[gvals]), "C1", label="EMCEE Fits")

    if args.symfit:
        ax[0].plot(wave[gvals], g21_fit_y, "g-", label="Sym Best Fit", alpha=0.7)
    ax[0].plot(wave[gvals], g21_asym_fit_y, "b-", label="Best Fit")

    mmy = np.array([min(g21_fit_y), max(g21_fit_y)])
    if obsext.type == "elx":
        mmy[0] = min(
            [mmy[0], -1.0 * g21_fit.Av_1.value, -1.0 * g21_asym_fit.Av_1.value]
        )
    mmd = 0.1 * (mmy[1] - mmy[0])
    ax[0].set_ylim(mmy + np.array([-1.0, 1.0]) * mmd)
    ax[0].set_xlim(1.0, 40.0)
    ax[0].set_xscale("log")
    if not args.notitle:
        ax[0].set_title(file)

    g21_comps = g21_asym_fit.copy()
    if obsext.type == "elx":
        g21_comps[0].sil1_amp = 0.0
    else:
        g21_comps.sil1_amp = 0.0
    ax[0].plot(wave[gvals], g21_comps(wave[gvals]), "k--", alpha=0.5)

    g21_comps = g21_asym_fit.copy()
    if obsext.type == "elx":
        g21_comps[0].sil2_amp = 0.0
    else:
        g21_comps.sil2_amp = 0.0
    ax[0].plot(wave[gvals], g21_comps(wave[gvals]), "k--", alpha=0.5)

    g21_comps = g21_asym_fit.copy()
    if obsext.type == "elx":
        g21_comps[0].sil1_amp = 0.0
        g21_comps[0].sil2_amp = 0.0
    else:
        g21_comps.sil1_amp = 0.0
        g21_comps.sil2_amp = 0.0
    ax[0].plot(wave[gvals], g21_comps(wave[gvals]), "k--", alpha=0.5)

    ax[0].legend(loc="best")

    # residuals
    ax[1].plot(wave[gvals], np.zeros((len(wave[gvals]))), "k--")
    if args.symfit:
        ax[1].plot(wave[gvals], y[gvals] - g21_fit_y, "g-", alpha=0.7)

    gbands = obsext.waves["BAND"] > (1.0 * u.micron)
    ax[1].errorbar(
        obsext.waves["BAND"][gbands].value,
        obsext.exts["BAND"][gbands] - g21_asym_fit(obsext.waves["BAND"][gbands]),
        yerr=obsext.uncs["BAND"][gbands],
        fmt="bo",
        mfc="white",
    )
    ax[1].plot(
        obsext.waves["IRS"].value,
        obsext.exts["IRS"] - g21_asym_fit(obsext.waves["IRS"]),
        "b-",
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
