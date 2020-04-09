import matplotlib.pyplot as plt
import numpy as np
import argparse
import warnings

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
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    # get a saved extnction curve
    file = args.extfile
    # file = '/home/kgordon/Python_git/spitzer_mir_ext/fits/hd147889_hd064802_ext.fits'
    ofile = file.replace(".fits", "_fm90.fits")
    ext = ExtData(filename=file)
    ext.trans_elv_alav(av=float(ext.columns["AV"][0]))
    gindxs = ext.npts["IUE"] > 0
    x = 1.0 / ext.waves["IUE"][gindxs].to(u.micron).value
    y = ext.exts["IUE"][gindxs]
    y_unc = ext.uncs["IUE"][gindxs]
    gindxs = (x > 3.5) & (x < 8.0)

    # initialize the model
    fm90_init = FM90()

    # pick the fitter
    fit = LevMarLSQFitter()
    # fit2 = SPyMinimizeFitter()
    nsteps = args.nsteps
    fit3 = EmceeFitter(nsteps=nsteps)

    # fit the data to the FM90 model using the fitter
    #   use the initialized model as the starting point
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        fm90_fit = fit(fm90_init, x[gindxs], y[gindxs], weights=1.0 / y_unc[gindxs])
        # fm90_fit2 = fit2(fm90_init, x[gindxs], y[gindxs], weights=1.0 / y_unc[gindxs])
        fm90_fit3 = fit3(fm90_fit, x[gindxs], y[gindxs], weights=1.0 / y_unc[gindxs])

    # checking the uncertainties
    print("Best Fit Parameters")
    print(fm90_fit3.parameters)
    print("Symmetric uncertainties")
    print(fm90_fit3.uncs)
    print("Plus uncertainties")
    print(fm90_fit3.uncs_plus)
    print("Minus uncertainties")
    print(fm90_fit3.uncs_minus)

    for i, pname in enumerate(fm90_fit3.param_names):
        # now set uncertainties on the parameter objects themselves
        param = getattr(fm90_fit3, pname)
        if param.posterior is not None:
            print(
                "posterior: ",
                pname,
                param.posterior.pdf_mean(),
                param.posterior.pdf_std(),
            )
        print("parameter: ", pname, param.value, param.unc)

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

    # print(fit3.fit_info['perparams'])

    # print(fit3.fit_info['sampler'].get_autocorr_time())

    # plot the observed data, initial guess, and final fit
    fig, ax = plt.subplots()

    # remove pesky x without units warnings
    x /= u.micron

    # ax.errorbar(x, y, yerr=y_unc[gindxs], fmt='ko', label='Observed Curve')
    # ax.plot(x[gindxs], fm90_init(x[gindxs]), label='Initial guess')
    ax.plot(x, y, label="Observed Curve")
    ax.plot(x[gindxs], fm90_fit3(x[gindxs]), label="emcee")
    # ax.plot(x[gindxs], fm90_fit2(x[gindxs]), label="scipy.minimize")
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
    # ax.set_ylabel('$A(x)/A(V)$')
    ax.set_ylabel(r"$A(\lambda)/A(V)$")

    ax.set_title(file)
    # ax.set_title('FM90 Fit to G09_MWAvg curve')

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
