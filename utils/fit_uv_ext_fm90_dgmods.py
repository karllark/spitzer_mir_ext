import matplotlib.pyplot as plt
import numpy as np
import argparse
import warnings

# from multiprocessing import Pool
# using a Pool does not work in this setup it seems

from astropy.modeling.fitting import LevMarLSQFitter
import astropy.units as u

from dust_extinction.shapes import FM90

from dust_extinction.grain_models import D03, ZDA04, J13


if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dgmodel",
        help="dust grain model to use",
        choices=["D03", "ZDA04", "J13"],
        default="D03",
    )
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    wave = np.logspace(np.log10(0.1), np.log10(0.4), num=100) * u.micron
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
    x = 1. / wave

    gindxs = (x > (3.3 / u.micron)) & (x < (8.0 / u.micron))
    # gindxs = ((x > 3.3) & (x < 8.1)) | ((x > 8.35) & (x < 8.6))

    # initialize the model
    fm90_init = FM90(C3=0.75)

    fm90_init.C1.bounds = (0., 3.)
    fm90_init.C2.bounds = (-0.1, 0.6)
    fm90_init.C3.bounds = (0., 2.5)
    fm90_init.C4.bounds = (0., 1.)
    fm90_init.xo.bounds = (4.5, 4.9)
    fm90_init.gamma.bounds = (0.6, 1.5)

    # Set up the backend to save the samples for the emcee runs
    emcee_samples_file = ofile.replace(".fits", ".h5")

    # pick the fitter
    fit = LevMarLSQFitter()

    # fit the data to the FM90 model using the fitter
    #   use the initialized model as the starting point
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        fm90_fit = fit(fm90_init, x[gindxs].value, y[gindxs])

    # setup the parameters for saving
    fm90_best_params = (fm90_fit.param_names, fm90_fit.parameters)
    print(fm90_best_params)

    # plot the observed data, initial guess, and final fit
    fig, ax = plt.subplots()

    # remove pesky x without units warnings
    x /= u.micron

    # ax.plot(x[gindxs], fm90_init(x[gindxs]), label='Initial guess')
    ax.plot(x, y, label="Observed Curve")
    ax.plot(x[gindxs], fm90_fit(x[gindxs].value), label="LevMarLSQ")

    ax.set_xlabel(r"$x$ [$\mu m^{-1}$]")
    ax.set_ylabel(r"$A(\lambda)/A(V)$")

    ax.set_title(args.dgmodel)

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
