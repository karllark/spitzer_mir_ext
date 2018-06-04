from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import argparse
import warnings
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

import emcee
import corner

from astropy.modeling.fitting import LevMarLSQFitter
from astropy.utils.exceptions import AstropyWarning

from dust_extinction.shapes import P92
from dust_extinction.conversions import AxAvToExv
from measure_extinction.stardata import StarData
from measure_extinction.extdata import ExtData


def tie_amps_SIL2_to_SIL1(model):
    """
    function to tie the SIL2 amplitude to the SIL1 amplitude
    """
    return model.SIL1_amp_0


def tie_amps_FIR_to_SIL1(model):
    """
    function to tie the FIR amplitude to the SIL1 amplitude
    """
    return (0.012/0.002)*model.SIL1_amp_0


class P92_Elv(P92 | AxAvToExv):
    """
    Evalute P92 on E(x-V) data including solving for A(V)
    """


def lnprob(params, x, y, uncs, model, param_names):
    """
    Log likelihood

    Parameters
    ----------
    params : array of floats
        parameters for evaluting the model

    x, y, uncs : array
        x, y, and y uncertainties for the data

    model : astropy model
        model for evaluation

    param_name : array of str
        names of the model parameters to update based on params
        should not include any fixed parameters
    """
    param_dict = dict(zip(model.param_names, model.parameters))
    for k, cname in enumerate(param_names):
        # impose the bounds
        if model.bounds[cname][0] is not None:
            if params[k] < model.bounds[cname][0]:
                return -np.inf
        if model.bounds[cname][1] is not None:
            if params[k] > model.bounds[cname][1]:
                return -np.inf
        # otherwise, set the requested value
        exec("model.{}.value = {}".format(cname, params[k]))
    return -0.5*np.sum(((model(x) - y)/uncs)**2)


def get_best_fit_params(sampler):
    """
    Get the best fit parameters from all the walkers

    Parameters
    ----------
    sample : emcee sampler object

    Returns
    -------
    fit_params_best : array of floats
        parameters of the best fit (largest likelihood sample)
    """
    max_lnp = -1e32
    nwalkers = len(sampler.lnprobability)
    fit_params_best = None
    for k in range(nwalkers):
        tmax_lnp = np.max(sampler.lnprobability[k])
        if tmax_lnp > max_lnp:
            max_lnp = tmax_lnp
            indxs, = np.where(sampler.lnprobability[k] == tmax_lnp)
            fit_params_best = sampler.chain[k, indxs[0], :]
    return fit_params_best


def p92_emcee(x, y, uncs,
              model, fit_param_names=None,
              threads=1, return_sampler=False,
              nburn=100, nsteps=500):
    """
    Fit the model using the emcee MCMC sampler

    Parameters
    ----------
    x, y, uncs : array
        x, y, and y uncertainties for the data

    model : astropy model
        model to fit
        starting position taken from model paramter values

    fit_param_names : list of string, optional
        list of parameters to fit
        default is to fit all non-fixed parameters

    threads : int
        number of threads to use for MCMC run

    nburn : int
        number of steps for the MCMC burn in

    nsteps : int
        number of steps for the MCMC sampling

    return_sampler: booelean
        return emcee sampler, return is now (best_fit_model, sampler)
    """

    model_copy = model.copy()

    # get a list of non-fixed parameters
    if fit_param_names is None:
        fit_param_names = []
        for cname in model_copy.param_names:
            if not model_copy.fixed[cname]:
                fit_param_names.append(cname)

    # sampler setup
    ndim = len(fit_param_names)
    nwalkers = 10*ndim

    # needed for priors
    # model_copy.bounds

    # inital guesses at parameters
    p0_list = []
    param_dict = dict(zip(model_copy.param_names, model_copy.parameters))
    for cname in fit_param_names:
        p0_list.append(param_dict[cname])
    p0 = np.array(p0_list)

    # check if any parameters are zero and make them sligthly larger
    p0[p0 == 0.0] = 2.4e-3

    print(fit_param_names)
    print(p0)

    # setting up the walkers to start "near" the inital guess
    p = [p0*(1+1e-4*np.random.normal(0, 1., ndim)) for k in range(nwalkers)]

    # for the parameters with min/max bounds set ("good priors")
    # sample from the prior
#    for k, cname in enumerate(fit_param_names):
#        if ((model.bounds[cname][0] is not None)
#                & (model.bounds[cname][1] is not None)):
#            svals = np.random.uniform(model.bounds[cname][0],
#                                      model.bounds[cname][1],
#                                      nwalkers)
#            for i in range(nwalkers):
#                p[i][k] = svals[i]
#            print(p)

    # ensure all the walkers start within the bounds
    param_dict = dict(zip(model.param_names, model.parameters))
    for cp in p:
        for k, cname in enumerate(fit_param_names):
            # check the bounds
            if model.bounds[cname][0] is not None:
                if cp[k] < model.bounds[cname][0]:
                    cp[k] = model.bounds[cname][0]
                    # print('min: ', cname, cp[k])
            if model.bounds[cname][1] is not None:
                if cp[k] > model.bounds[cname][1]:
                    cp[k] = model.bounds[cname][1]
                    # print('max: ', cname, cp[k])

    # setup the sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, threads=threads,
                                    args=(x, y, uncs, model_copy,
                                          fit_param_names))

    if nburn is not None:
        # burn in the walkers
        pos, prob, state = sampler.run_mcmc(p, nburn)
        # rest the sampler
        sampler.reset()

    # do the full sampling
    pos, prob, state = sampler.run_mcmc(pos, nsteps, rstate0=state)

    # best fit parameters
    best_params = get_best_fit_params(sampler)

    # percentile parameters
    samples = sampler.chain.reshape((-1, ndim))
    per_params = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                     zip(*np.percentile(samples, [16, 50, 84],
                                        axis=0)))
    # set the returned model parameters to the p50 values
    for k, val in enumerate(per_params):
        exec("model_copy.{}.value = {}".format(fit_param_names[k],
                                               val[0]))
        print(fit_param_names[k], best_params[k], val)

    clean_pnames = [pname[:-2] for pname in fit_param_names]
    model_copy.p92_emcee_param_names = clean_pnames
    model_copy.p92_emcee_per_params = per_params

    if return_sampler:
        return (model_copy, sampler)
    else:
        return model_copy


def plot_emcee_results(sampler,
                       fit_param_names,
                       filebase=''):
    """
    Plot the standard triangle and diagnostic walker plots
    """

    # plot the walker chains for all parameters
    nwalkers, nsteps, ndim = sampler.chain.shape
    fig, ax = plt.subplots(ndim, sharex=True, figsize=(13, 13))
    walk_val = np.arange(nsteps)
    for i in range(ndim):
        for k in range(nwalkers):
            ax[i].plot(walk_val, sampler.chain[k, :, i], '-')
            ax[i].set_ylabel(fit_param_names[i])
    fig.savefig('%s_walker_param_values.png' % filebase)
    plt.close(fig)

    # plot the 1D and 2D likelihood functions in a traditional triangle plot
    samples = sampler.chain.reshape((-1, ndim))
    fig = corner.corner(samples, labels=fit_param_names, show_titles=True,
                        title_fmt='.3f', use_math_text=True)
    fig.savefig('%s_param_triangle.png' % filebase)
    plt.close(fig)


if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("redstarname", help="name of reddened star")
    parser.add_argument("compstarname", help="name of comparision star")
    parser.add_argument("--path", help="base path to observed data",
                        default="/home/kgordon/Dust/Ext/")
    parser.add_argument("--emcee", help="run EMCEE fit",
                        action="store_true")
    parser.add_argument("--nburn", type=int, default=100,
                        help='# of burn steps in MCMC chain')
    parser.add_argument("--nsteps", type=int, default=500,
                        help='# of steps in MCMC chain')
    parser.add_argument("--threads", type=int, default=1,
                        help="number of threads for EMCEE run")
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

    # output filebase
    filebase = 'fits/%s_%s' % (args.redstarname, args.compstarname)

    # calculate the extinction curve
    extdata = ExtData()
    extdata.calc_elv(redstarobs, compstarobs)

    # fit the calculated extinction curve
    # get an observed extinction curve to fit
    (x, y, y_unc) = extdata.get_fitdata(['BAND', 'IUE', 'IRS'],
                                        remove_uvwind_region=True,
                                        remove_lya_region=True)

    # determine the initial guess at the A(V) values
    #  just use the average at wavelengths > 5
    #  limit as lambda -> inf, E(lamda-V) -> -A(V)
    indxs, = np.where(1./x > 5.0)
    av_guess = -1.0*np.average(y[indxs])
    if not np.isfinite(av_guess):
        av_guess = 1.0

    # initialize the model
    #    a few tweaks to the starting parameters helps find the solution
    p92_init = P92_Elv(BKG_amp_0=200.,
                       FUV_amp_0=100.0, FUV_lambda_0=0.06,
                       # FIR_amp_0 = 0.0,
                       Av_1=av_guess)

    # fix a number of the parameters
    #   mainly to avoid fitting parameters that are constrained at
    #   wavelengths where the observed data for this case does not exist
    p92_init.BKG_b_0.fixed = True
    p92_init.FUV_lambda_0.fixed = True
    p92_init.FUV_b_0.fixed = True
    p92_init.FUV_n_0.fixed = True
    p92_init.NUV_b_0.fixed = True
    p92_init.SIL1_b_0.fixed = True
    p92_init.SIL2_lambda_0.fixed = True
    p92_init.SIL2_b_0.fixed = True
    # p92_init.FIR_amp_0.fixed = True
    p92_init.FIR_lambda_0.fixed = True
    p92_init.FIR_b_0.fixed = True

    p92_init.Av_1.bounds = [0.1, None]

    # p92_init.SIL2_amp_0.tied = tie_amps_SIL2_to_SIL1
    # p92_init.FIR_amp_0.tied = tie_amps_FIR_to_SIL1

    # pick the fitter
    fit = LevMarLSQFitter()

    # fit the data to the P92 model using the fitter
    # p92_fit = fit(p92_init, x, y)
    p92_fit = fit(p92_init, x, y, weights=1.0/y_unc)

    clean_pnames = [pname[:-2] for pname in p92_init.param_names]
    print("P92 names")
    print(clean_pnames)
    print("initital parameters")
    print(p92_init.parameters)
    print("best fit parameters")
    print(p92_fit.parameters)

    best_fit_Av = p92_fit.Av_1.value

    if args.emcee:
        # run the emcee fitter to get proper fit parameter uncertainties
        fit_param_names = ['Av_1',
                           'BKG_amp_0', 'BKG_lambda_0',
                           'NUV_amp_0', 'NUV_lambda_0',
                           'FUV_amp_0',
                           'SIL1_amp_0', 'SIL1_lambda_0',
                           'SIL2_amp_0',
                           'FIR_amp_0']
        emcee_results = p92_emcee(x, y, y_unc, p92_fit,
                                  nburn=args.nburn,
                                  nsteps=args.nsteps,
                                  threads=args.threads,
                                  fit_param_names=fit_param_names,
                                  return_sampler=True)
        p92_fit_emcee, sampler = emcee_results
        clean_pnames_emcee = [pname[:-2] for pname in fit_param_names]
        plot_emcee_results(sampler, clean_pnames_emcee, filebase=filebase)

    # save the extinction curve and fit
    warnings.simplefilter('ignore', category=AstropyWarning)
    out_fname = "fits/%s_%s_ext.fits" % (args.redstarname, args.compstarname)
    extdata.save_ext_data(out_fname,
                          p92_best_params=(clean_pnames, p92_fit.parameters))

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
    fig, ax = plt.subplots(figsize=(12, 8))

    # subplot
    ax2 = plt.axes([.60, .35, .35, .35])

    # plot the bands and all spectra for this star
    extdata.plot_ext(ax, color='k', alpha=0.5)
    extdata.plot_ext(ax2, color='k', alpha=0.5)

    # ax.plot(1./x, p92_init(x), 'p-', label='Initial guess', alpha=0.25)
    # ax2.plot(1./x, p92_init(x), 'p-', label='Initial guess', alpha=0.25)
    if args.emcee:
        # plot the walkers with transparancy
        nwalkers, nsteps, ndim = sampler.chain.shape
        nwalkers = 1
        # nsteps = 100
        for k in range(nwalkers):
            for j in range(nsteps):
                for i, pname in enumerate(fit_param_names):
                    exec("p92_fit_emcee.{}.value = {}".format(pname,
                                                sampler.chain[k, j, i]))
                ax.plot(1./x, p92_fit_emcee(x), 'b-', alpha=0.01)
                ax2.plot(1./x, p92_fit_emcee(x), 'b-', alpha=0.01)

        # p50 results
        ax.plot(1./x, p92_fit_emcee(x), 'b-', label='EMCEE Fits')
        ax2.plot(1./x, p92_fit_emcee(x), 'b-')
    ax.plot(1./x, p92_fit(x), 'r-', label='Best Fit')
    ax2.plot(1./x, p92_fit(x), 'r-')

    # show components of best fitter
    p92_comps = p92_fit.copy()
    p92_comps.FUV_amp_0 = 0.0
    p92_comps.NUV_amp_0 = 0.0
    p92_comps.SIL1_amp_0 = 0.0
    p92_comps.SIL2_amp_0 = 0.0
    p92_comps.FIR_amp_0 = 0.0
    ax.plot(1./x, p92_comps(x), 'k-', alpha=0.5)
    ax2.plot(1./x, p92_comps(x), 'k-', alpha=0.5)

    ax.plot(1./x, best_fit_Av*np.full((len(x)), -1.0), '-', label='-A(V)')
    ax2.plot(1./x, best_fit_Av*np.full((len(x)), -1.0), '-')

    p92_comps = p92_fit.copy()
    p92_comps.NUV_amp_0 = 0.0
    p92_comps.SIL1_amp_0 = 0.0
    p92_comps.SIL2_amp_0 = 0.0
    p92_comps.FIR_amp_0 = 0.0
    ax.plot(1./x, p92_comps(x), 'k--', alpha=0.5)
    ax2.plot(1./x, p92_comps(x), 'k--', alpha=0.5)

    p92_comps = p92_fit.copy()
    p92_comps.FUV_amp_0 = 0.0
    p92_comps.SIL1_amp_0 = 0.0
    p92_comps.SIL2_amp_0 = 0.0
    p92_comps.FIR_amp_0 = 0.0
    ax.plot(1./x, p92_comps(x), 'k--', alpha=0.5)
    ax2.plot(1./x, p92_comps(x), 'k--', alpha=0.5)

    p92_comps = p92_fit.copy()
    p92_comps.FUV_amp_0 = 0.0
    p92_comps.NUV_amp_0 = 0.0
    p92_comps.SIL2_amp_0 = 0.0
    p92_comps.FIR_amp_0 = 0.0
    ax.plot(1./x, p92_comps(x), 'k--', alpha=0.5)
    ax2.plot(1./x, p92_comps(x), 'k--', alpha=0.5)

    p92_comps = p92_fit.copy()
    p92_comps.FUV_amp_0 = 0.0
    p92_comps.NUV_amp_0 = 0.0
    p92_comps.SIL1_amp_0 = 0.0
    p92_comps.FIR_amp_0 = 0.0
    ax.plot(1./x, p92_comps(x), 'k--', alpha=0.5)
    ax2.plot(1./x, p92_comps(x), 'k--', alpha=0.5)

    p92_comps = p92_fit.copy()
    p92_comps.FUV_amp_0 = 0.0
    p92_comps.NUV_amp_0 = 0.0
    p92_comps.SIL1_amp_0 = 0.0
    p92_comps.SIL2_amp_0 = 0.0
    ax.plot(1./x, p92_comps(x), 'k--', alpha=0.5)
    ax2.plot(1./x, p92_comps(x), 'k--', alpha=0.5)

    # finish configuring the plot
    ax.set_yscale('linear')
    ax.set_xscale('log')
    ax.set_xlabel('$\lambda$ [$\mu m$]', fontsize=1.3*fontsize)
    ax.set_ylabel(extdata._get_ext_ytitle(extdata.type),
                  fontsize=1.3*fontsize)
    ax.tick_params('both', length=10, width=2, which='major')
    ax.tick_params('both', length=5, width=1, which='minor')
    ax.legend()

    # finish configuring the subplot
    sp_xlim = [2.0, 35.0]
    ax2.set_xlim(sp_xlim)
    # ax2.set_ylim(-best_fit_Av-0.1, -best_fit_Av+0.5)
    indxs, = np.where((x > 1.0/sp_xlim[1]) & (x < 1.0/sp_xlim[0]))
    ax2.set_ylim(min([min(p92_fit(x)[indxs]), -best_fit_Av])-0.1,
                 max(p92_fit(x)[indxs])+0.1)

    # use the whitespace better
    warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")
    fig.tight_layout()

    # plot or save to a file
    outname = "%s_ext" % filebase
    if args.png:
        fig.savefig(outname+'.png')
    elif args.eps:
        fig.savefig(outname+'.eps')
    elif args.pdf:
        fig.savefig(outname+'.pdf')
    else:
        plt.show()
