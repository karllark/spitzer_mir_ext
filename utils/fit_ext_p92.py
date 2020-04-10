import argparse
import warnings
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

import astropy.units as u
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.modeling.fitting import _fitter_to_model_params
from astropy.utils.exceptions import AstropyWarning

# from dust_extinction.shapes import P92
from dust_extinction.conversions import AxAvToExv
from measure_extinction.extdata import ExtData

from models_mcmc_extension import EmceeFitter
from P92_mod import P92_mod


def tie_amps_SIL2_to_SIL1(model):
    """
    function to tie the SIL2 amplitude to the SIL1 amplitude
    """
    return model.SIL1_amp_0


def tie_amps_FIR_to_SIL1(model):
    """
    function to tie the FIR amplitude to the SIL1 amplitude
    """
    return (0.012 / 0.002) * model.SIL1_amp_0


def clean_pnames(pnames):
    """
    function to clean of the _? part of the names due to making a CompoundModel
    """
    clean_pnames = [cpname[:-2] for cpname in pnames]
    return clean_pnames


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
    parser.add_argument(
        "--nthreads", type=int, default=1, help="number of threads for MCMC run"
    )
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    # get a saved extnction curve
    file = args.extfile
    # file = '/home/kgordon/Python_git/spitzer_mir_ext/fits/hd147889_hd064802_ext.fits'
    ofile = file.replace(".fits", "_P92.fits")
    extdata = ExtData(filename=file)

    # get an observed extinction curve to fit
    (wave, y, y_unc) = extdata.get_fitdata(
        ["BAND", "IUE", "IRS"], remove_uvwind_region=True, remove_lya_region=True
    )
    # ["BAND", "IUE", "IRS"], remove_uvwind_region=True, remove_lya_region=True
    # remove data affected by Ly-alpha absorption/emission
    gindxs = wave > (1.0 / 8.0) * u.micron
    wave = wave[gindxs]
    y = y[gindxs]
    y_unc = y_unc[gindxs]

    # remove units as fitting routines often cannot take numbers with units
    x = wave.to(1.0 / u.micron, equivalencies=u.spectral()).value

    # determine the initial guess at the A(V) values
    #  just use the average at wavelengths > 5
    #  limit as lambda -> inf, E(lamda-V) -> -A(V)
    (indxs,) = np.where(1.0 / x > 5.0)
    av_guess = -1.0 * np.average(y[indxs])
    if not np.isfinite(av_guess):
        av_guess = 1.0

    # initialize the model
    #    a few tweaks to the starting parameters helps find the solution
    p92_init = P92_mod(BKG_amp=200.0, FUV_amp=100.0, FUV_lambda=0.06) | AxAvToExv(
        Av=av_guess
    )

    # fix a number of the parameters
    #   mainly to avoid fitting parameters that are constrained at
    #   wavelengths where the observed data for this case does not exist
    p92_init.BKG_lambda_0.fixed = True
    p92_init.BKG_width_0.fixed = True
    # p92_init.FUV_amp_0.fixed = True
    p92_init.FUV_lambda_0.fixed = True
    p92_init.FUV_b_0.fixed = True
    p92_init.FUV_n_0.fixed = True
    # p92_init.NUV_amp_0.fixed = True
    # p92_init.NUV_lambda_0.fixed = True
    # p92_init.NUV_width_0.fixed = True
    p92_init.SIL2_lambda_0.fixed = True
    p92_init.SIL2_width_0.fixed = True
    p92_init.FIR_lambda_0.fixed = True
    p92_init.FIR_width_0.fixed = True

    p92_init.Av_1.bounds = [0.1, None]

    # extra for HD 147933 (4 Nov 2019)
    # p92_init.SIL1_amp_0.bounds = [0.005, None]
    # p92_init.BKG_lambda_0.fixed = False
    # p92_init.BKG_width_0.fixed = False
    # p92_init.FUV_lambda_0.fixed = False

    # p92_init.SIL2_amp_0.tied = tie_amps_SIL2_to_SIL1
    # p92_init.FIR_amp_0.tied = tie_amps_FIR_to_SIL1

    # pick the fitter
    fit = LevMarLSQFitter()
    nsteps = args.nsteps
    fit2 = EmceeFitter(nsteps=nsteps, burnfrac=args.burnfrac)

    # fit the data to the P92 model using the fitter
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        p92_fit = fit(p92_init, x, y, weights=1.0 / y_unc, maxiter=10000, epsilon=0.001)
        p92_fit2 = fit2(p92_fit, x, y, weights=1.0 / y_unc)

    p92_best_params = (clean_pnames(p92_fit.param_names), p92_fit.parameters)
    p92_per_param_vals = zip(
        p92_fit2.parameters, p92_fit2.uncs_plus, p92_fit2.uncs_minus
    )
    p92_per_params = (clean_pnames(p92_fit2.param_names), list(p92_per_param_vals))

    # save the extinction curve and fit
    warnings.simplefilter("ignore", category=AstropyWarning)
    extdata.save(
        ofile, p92_best_params=p92_best_params, p92_per_params=p92_per_params
    )

    # make the standard mcmc plots
    fit2.plot_emcee_results(p92_fit2, filebase=ofile.replace(".fits", ""))

    best_fit_Av = p92_fit.Av_1.value
    x /= u.micron

    # plotting setup for easier to read plots
    # fontsize = 18
    fontsize = 10
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
    # fig, ax = plt.subplots(figsize=(6, 4))

    # subplot
    ax2 = plt.axes([0.60, 0.35, 0.35, 0.35])

    # plot the bands and all spectra for this star
    extdata.plot(ax, color="k", alpha=0.5)
    extdata.plot(ax2, color="k", alpha=0.5)

    # plot samples from the mcmc chaing
    flat_samples = fit2.fit_info["sampler"].get_chain(
        discard=int(0.1 * nsteps), flat=True
    )
    inds = np.random.randint(len(flat_samples), size=100)
    model_copy = p92_fit2.copy()
    for ind in inds:
        sample = flat_samples[ind]
        _fitter_to_model_params(model_copy, sample)
        ax.plot(1.0 / x, model_copy(x), "C1", alpha=0.05)
        ax2.plot(1.0 / x, model_copy(x), "C1", alpha=0.05)
    # for the figure legend
    ax.plot(
        1.0 / x, model_copy(x), "C1", alpha=0.05, label="EMCEE Fits"
    )
    ax2.plot(
        1.0 / x, model_copy(x), "C1", alpha=0.05, label="EMCEE Fits"
    )

    # ax.plot(1.0 / x, p92_init(x), "r--", label="P92 Init")
    ax.plot(1.0 / x, p92_fit(x), "r-", label="P92 Best Fit")
    ax2.plot(1.0 / x, p92_fit(x), "r-")

    # show components of best fitter
    p92_comps = p92_fit.copy()
    p92_comps.FUV_amp_0 = 0.0
    p92_comps.NUV_amp_0 = 0.0
    p92_comps.SIL1_amp_0 = 0.0
    p92_comps.SIL2_amp_0 = 0.0
    p92_comps.FIR_amp_0 = 0.0
    ax.plot(1.0 / x, p92_comps(x), "k-", alpha=0.5)
    ax2.plot(1.0 / x, p92_comps(x), "k-", alpha=0.5)

    ax.plot(1.0 / x, best_fit_Av * np.full((len(x)), -1.0), "-", label="-A(V)")
    ax2.plot(1.0 / x, best_fit_Av * np.full((len(x)), -1.0), "-")

    p92_comps = p92_fit.copy()
    p92_comps.NUV_amp_0 = 0.0
    p92_comps.SIL1_amp_0 = 0.0
    p92_comps.SIL2_amp_0 = 0.0
    p92_comps.FIR_amp_0 = 0.0
    ax.plot(1.0 / x, p92_comps(x), "k--", alpha=0.5)
    ax2.plot(1.0 / x, p92_comps(x), "k--", alpha=0.5)

    p92_comps = p92_fit.copy()
    p92_comps.FUV_amp_0 = 0.0
    p92_comps.SIL1_amp_0 = 0.0
    p92_comps.SIL2_amp_0 = 0.0
    p92_comps.FIR_amp_0 = 0.0
    ax.plot(1.0 / x, p92_comps(x), "k--", alpha=0.5)
    ax2.plot(1.0 / x, p92_comps(x), "k--", alpha=0.5)

    p92_comps = p92_fit.copy()
    p92_comps.FUV_amp_0 = 0.0
    p92_comps.NUV_amp_0 = 0.0
    p92_comps.SIL2_amp_0 = 0.0
    p92_comps.FIR_amp_0 = 0.0
    ax.plot(1.0 / x, p92_comps(x), "k--", alpha=0.5)
    ax2.plot(1.0 / x, p92_comps(x), "k--", alpha=0.5)

    p92_comps = p92_fit.copy()
    p92_comps.FUV_amp_0 = 0.0
    p92_comps.NUV_amp_0 = 0.0
    p92_comps.SIL1_amp_0 = 0.0
    p92_comps.FIR_amp_0 = 0.0
    ax.plot(1.0 / x, p92_comps(x), "k--", alpha=0.5)
    ax2.plot(1.0 / x, p92_comps(x), "k--", alpha=0.5)

    p92_comps = p92_fit.copy()
    p92_comps.FUV_amp_0 = 0.0
    p92_comps.NUV_amp_0 = 0.0
    p92_comps.SIL1_amp_0 = 0.0
    p92_comps.SIL2_amp_0 = 0.0
    ax.plot(1.0 / x, p92_comps(x), "k--", alpha=0.5)
    ax2.plot(1.0 / x, p92_comps(x), "k--", alpha=0.5)

    # finish configuring the plot
    ax.set_yscale("linear")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\lambda$ [$\mu m$]", fontsize=1.3 * fontsize)
    ax.set_ylabel(extdata._get_ext_ytitle(extdata.type), fontsize=1.3 * fontsize)
    ax.tick_params("both", length=10, width=2, which="major")
    ax.tick_params("both", length=5, width=1, which="minor")
    ax.legend()

    # finish configuring the subplot
    sp_xlim = [2.0, 35.0]
    ax2.set_xlim(sp_xlim)
    # ax2.set_ylim(-best_fit_Av-0.1, -best_fit_Av+0.5)
    (indxs,) = np.where((x.value > 1.0 / sp_xlim[1]) & (x.value < 1.0 / sp_xlim[0]))
    ax2.set_ylim(
        min([min(p92_fit(x)[indxs]), -best_fit_Av]) - 0.1, max(p92_fit(x)[indxs]) + 0.1
    )

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
