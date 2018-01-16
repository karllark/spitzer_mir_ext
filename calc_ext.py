from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import argparse
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

from astropy.modeling.fitting import LevMarLSQFitter

from dust_extinction.dust_extinction import (P92, AxAvToExv)
from measure_extinction.stardata import StarData
from measure_extinction.extdata import ExtData


class P92_Elv(P92 | AxAvToExv):
    """
    Evalute P92 on E(x-V) data including solving for A(V)
    """


if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("redstarname", help="name of reddened star")
    parser.add_argument("compstarname", help="name of comparision star")
    parser.add_argument("--path", help="base path to observed data",
                        default="/home/kgordon/Dust/Ext/")
    parser.add_argument("--png", help="save figure as a png file",
                        action="store_true")
    parser.add_argument("--eps", help="save figure as an eps file",
                        action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file",
                        action="store_true")
    args = parser.parse_args()

    # read in the observed data for both stars
    redstarobs = StarData('DAT_files/%s.dat' % args.redstarname,
                          path=args.path)
    compstarobs = StarData('DAT_files/%s.dat' % args.compstarname,
                           path=args.path)

    # calculate the extinction curve
    extdata = ExtData()
    extdata.calc_elv(redstarobs, compstarobs)

    # fit the calculated extinction curve
    # get an observed extinction curve to fit
    iue_gindxs = np.where(extdata.uncs['IUE'] > 0.)
    irs_gindxs = np.where(extdata.uncs['IRS'] > 0.)
    x = np.concatenate([extdata.waves['IUE'][iue_gindxs],
                        extdata.waves['IRS'][irs_gindxs],
                        extdata.waves['BAND']])
    y = np.concatenate([extdata.exts['IUE'][iue_gindxs],
                        extdata.exts['IRS'][irs_gindxs],
                        extdata.exts['BAND']])
    y_unc = np.concatenate([extdata.uncs['IUE'][iue_gindxs],
                            extdata.uncs['IRS'][irs_gindxs],
                            extdata.uncs['BAND']])

    # sort the data
    sindxs = np.argsort(x)
    x = 1.0/x[sindxs]
    y = y[sindxs]
    y_unc = y_unc[sindxs]

    # determine the initial guess at the A(V) values
    #  just use the average at wavelengths > 5
    #  limit as lambda -> inf, E(lamda-V) -> -A(V)
    indxs, = np.where(1./x > 5.0)
    av_guess = -1.0*np.average(y[indxs])

    # initialize the model
    #    a few tweaks to the starting parameters helps find the solution
    p92_init = P92_Elv(BKG_amp_0=200.,
                       FUV_amp_0=100.0, FUV_lambda_0=0.06,
                       Av_1=av_guess)

    # fix a number of the parameters
    #   mainly to avoid fitting parameters that are constrained at
    #   wavelengths where the observed data for this case does not exist
    p92_init.FUV_lambda_0.fixed = True
    # p92_init.FIR_amp_0.fixed = True
    p92_init.SIL2_lambda_0.fixed = True
    p92_init.FIR_lambda_0.fixed = True
    # p92_init.FIR_b_0.fixed = True

    # pick the fitter
    fit = LevMarLSQFitter()

    # fit the data to the P92 model using the fitter
    p92_fit = fit(p92_init, x, y, weights=1.0/y_unc)

    # print(fit.fit_info)
    print(p92_init._parameters)
    print(p92_fit._parameters)

    # plotting setup for easier to read plots
    fontsize = 18
    font = {'size': fontsize}
    matplotlib.rc('font', **font)
    matplotlib.rc('lines', linewidth=1)
    matplotlib.rc('axes', linewidth=2)
    matplotlib.rc('xtick.major', width=2)
    matplotlib.rc('xtick.minor', width=2)
    matplotlib.rc('ytick.major', width=2)
    matplotlib.rc('ytick.minor', width=2)

    # setup the plot
    fig, ax = plt.subplots(figsize=(10, 13))

    # plot the bands and all spectra for this star
    extdata.plot_ext(ax)

    ax.plot(1./x, p92_init(x), label='Initial guess')
    ax.plot(1./x, p92_fit(x), label='Fitted model')

    # finish configuring the plot
    ax.set_yscale('linear')
    ax.set_xscale('log')
    ax.set_xlabel('$\lambda$ [$\mu m$]', fontsize=1.3*fontsize)
    ax.set_ylabel(extdata._get_ext_ytitle(extdata.type),
                  fontsize=1.3*fontsize)
    ax.tick_params('both', length=10, width=2, which='major')
    ax.tick_params('both', length=5, width=1, which='minor')

    # use the whitespace better
    fig.tight_layout()

    # plot or save to a file
    save_str = '_spec'
    if args.png:
        fig.savefig(args.starname.replace('.dat', save_str+'.png'))
    elif args.eps:
        fig.savefig(args.starname.replace('.dat', save_str+'.eps'))
    elif args.pdf:
        fig.savefig(args.starname.replace('.dat', save_str+'.pdf'))
    else:
        plt.show()
