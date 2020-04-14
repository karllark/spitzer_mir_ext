import matplotlib.pyplot as plt
import numpy as np
import argparse
import warnings

# from multiprocessing import Pool
# using a Pool does not work in this setup it seems

from astropy.modeling.fitting import LevMarLSQFitter
from astropy.modeling.fitting import _fitter_to_model_params
import astropy.units as u

from models_mcmc_extension import EmceeFitter

from dust_extinction.shapes import FM90

from measure_extinction.extdata import ExtData


if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("extfile", help="file with extinction curve")
    parser.add_argument(
        "--nsteps", type=int, default=1000, help="# of steps in MCMC chain"
    )
    parser.add_argument(
        "--burnfrac", type=float, default=0.1, help="fraction of MCMC chain to burn"
    )
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    # get a saved extnction curve
    file = args.extfile
    # file = '/home/kgordon/Python_git/spitzer_mir_ext/fits/hd147889_hd064802_ext.fits'
    ofile = file.replace(".fits", "_FM90.fits")
    ext = ExtData(filename=file)
    ext.trans_elv_alav(av=float(ext.columns["AV"][0]))
    gindxs = ext.npts["IUE"] > 0
    x = 1.0 / ext.waves["IUE"][gindxs].to(u.micron).value
    y = ext.exts["IUE"][gindxs]
    y_unc = ext.uncs["IUE"][gindxs]
    gindxs = (x > 3.3) & (x < 8.0)
    # gindxs = ((x > 3.3) & (x < 8.1)) | ((x > 8.35) & (x < 8.6))

    # initialize the model
    fm90_init = FM90()

    # Set up the backend to save the samples for the emcee runs
    emcee_samples_file = ofile.replace(".fits", ".h5")

    # pick the fitter
    fit = LevMarLSQFitter()
    nsteps = args.nsteps
    fit3 = EmceeFitter(
        nsteps=nsteps, burnfrac=args.burnfrac, save_samples=emcee_samples_file
    )

    # modify weights to make sure the 2175 A bump is fit
    weights = 1.0 / y_unc[gindxs]
    # weights[(x[gindxs] < 5.9)] *= 4.0
    weights[(x[gindxs] > 4.0) & (x[gindxs] < 5.1)] *= 10.0

    # fit the data to the FM90 model using the fitter
    #   use the initialized model as the starting point
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        fm90_fit = fit(fm90_init, x[gindxs], y[gindxs], weights=weights)
        fm90_fit3 = fit3(fm90_fit, x[gindxs], y[gindxs], weights=weights)

    print("autocorr tau = ", fit3.fit_info["sampler"].get_autocorr_time(quiet=True))

    # setup the parameters for saving
    fm90_best_params = (fm90_fit.param_names, fm90_fit.parameters)
    fm90_per_param_vals = zip(
        fm90_fit3.parameters, fm90_fit3.uncs_plus, fm90_fit3.uncs_minus
    )
    fm90_per_params = (fm90_fit3.param_names, list(fm90_per_param_vals))

    # save extinction and fit parameters
    ext.save(ofile, fm90_best_params=fm90_best_params, fm90_per_params=fm90_per_params)

    # make the standard mcmc plots
    fit3.plot_emcee_results(fm90_fit3, filebase=ofile.replace(".fits", ""))

    # plot the observed data, initial guess, and final fit
    fig, ax = plt.subplots()

    # remove pesky x without units warnings
    x /= u.micron

    # ax.plot(x[gindxs], fm90_init(x[gindxs]), label='Initial guess')
    ax.plot(x, y, label="Observed Curve")
    ax.plot(x[gindxs], fm90_fit3(x[gindxs]), label="emcee")
    ax.plot(x[gindxs], fm90_fit(x[gindxs]), label="LevMarLSQ")

    # plot samples from the mcmc chaing
    flat_samples = fit3.fit_info["sampler"].get_chain(
        discard=int(0.1 * nsteps), flat=True
    )
    inds = np.random.randint(len(flat_samples), size=100)
    model_copy = fm90_fit3.copy()
    for ind in inds:
        sample = flat_samples[ind]
        _fitter_to_model_params(model_copy, sample)
        plt.plot(x[gindxs], model_copy(x[gindxs]), "C1", alpha=0.05)

    ax.set_xlabel(r"$x$ [$\mu m^{-1}$]")
    ax.set_ylabel(r"$A(\lambda)/A(V)$")

    ax.set_title(file)

    ax.legend(loc="best")
    plt.tight_layout()

    # plot or save to a file
    outname = ofile.replace(".fits", "")
    if args.png:
        fig.savefig(outname + ".png")
    elif args.pdf:
        fig.savefig(outname + ".pdf")
    else:
        plt.show()
