import argparse
import warnings
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

import emcee
import corner

import astropy.units as u
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.utils.exceptions import AstropyWarning
from astropy.modeling import Fittable1DModel, Parameter

from dust_extinction.shapes import P92
from dust_extinction.conversions import AxAvToExv
from dust_extinction.helpers import _get_x_in_wavenumbers, _test_valid_x_range
from measure_extinction.stardata import StarData
from measure_extinction.extdata import ExtData


class P92_mod(Fittable1DModel):
    """
    P92 extinction model calculation

    Parameters
    ----------
    BKG_amp : float
      background term amplitude
    BKG_lambda : float
      background term central wavelength
    BKG_b : float
      background term b coefficient
    BKG_n : float
      background term n coefficient [FIXED at n = 2]

    FUV_amp : float
      far-ultraviolet term amplitude
    FUV_lambda : float
      far-ultraviolet term central wavelength
    FUV_b : float
      far-ultraviolet term b coefficent
    FUV_n : float
      far-ultraviolet term n coefficient

    NUV_amp : float
      near-ultraviolet (2175 A) term amplitude
    NUV_lambda : float
      near-ultraviolet (2175 A) term central wavelength
    NUV_b : float
      near-ultraviolet (2175 A) term b coefficent
    NUV_n : float
      near-ultraviolet (2175 A) term n coefficient [FIXED at n = 2]

    SIL1_amp : float
      1st silicate feature (~10 micron) term amplitude
    SIL1_lambda : float
      1st silicate feature (~10 micron) term central wavelength
    SIL1_b : float
      1st silicate feature (~10 micron) term b coefficent
    SIL1_n : float
      1st silicate feature (~10 micron) term n coefficient [FIXED at n = 2]

    SIL2_amp : float
      2nd silicate feature (~18 micron) term amplitude
    SIL2_lambda : float
      2nd silicate feature (~18 micron) term central wavelength
    SIL2_b : float
      2nd silicate feature (~18 micron) term b coefficient
    SIL2_n : float
      2nd silicate feature (~18 micron) term n coefficient [FIXED at n = 2]

    FIR_amp : float
      far-infrared term amplitude
    FIR_lambda : float
      far-infrared term central wavelength
    FIR_b : float
      far-infrared term b coefficent
    FIR_n : float
      far-infrared term n coefficient [FIXED at n = 2]

    Notes
    -----
    P92 extinction model

    From Pei (1992)

    Applicable from the extreme UV to far-IR

    Example showing a P92 curve with components identified.

    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt
        import astropy.units as u

        from dust_extinction.shapes import P92

        fig, ax = plt.subplots()

        # generate the curves and plot them
        lam = np.logspace(-3.0, 3.0, num=1000)
        x = (1.0/lam)/u.micron

        ext_model = P92()
        ax.plot(1/x,ext_model(x),label='total')

        ext_model = P92(FUV_amp=0., NUV_amp=0.0,
                        SIL1_amp=0.0, SIL2_amp=0.0, FIR_amp=0.0)
        ax.plot(1./x,ext_model(x),label='BKG only')

        ext_model = P92(NUV_amp=0.0,
                        SIL1_amp=0.0, SIL2_amp=0.0, FIR_amp=0.0)
        ax.plot(1./x,ext_model(x),label='BKG+FUV only')

        ext_model = P92(FUV_amp=0.,
                        SIL1_amp=0.0, SIL2_amp=0.0, FIR_amp=0.0)
        ax.plot(1./x,ext_model(x),label='BKG+NUV only')

        ext_model = P92(FUV_amp=0., NUV_amp=0.0,
                        SIL2_amp=0.0)
        ax.plot(1./x,ext_model(x),label='BKG+FIR+SIL1 only')

        ext_model = P92(FUV_amp=0., NUV_amp=0.0,
                        SIL1_amp=0.0)
        ax.plot(1./x,ext_model(x),label='BKG+FIR+SIL2 only')

        ext_model = P92(FUV_amp=0., NUV_amp=0.0,
                        SIL1_amp=0.0, SIL2_amp=0.0)
        ax.plot(1./x,ext_model(x),label='BKG+FIR only')

        # Milky Way observed extinction as tabulated by Pei (1992)
        MW_x = [0.21, 0.29, 0.45, 0.61, 0.80, 1.11, 1.43, 1.82,
                2.27, 2.50, 2.91, 3.65, 4.00, 4.17, 4.35, 4.57, 4.76,
                5.00, 5.26, 5.56, 5.88, 6.25, 6.71, 7.18, 7.60,
                8.00, 8.50, 9.00, 9.50, 10.00]
        MW_x = np.array(MW_x)
        MW_exvebv = [-3.02, -2.91, -2.76, -2.58, -2.23, -1.60, -0.78, 0.00,
                     1.00, 1.30, 1.80, 3.10, 4.19, 4.90, 5.77, 6.57, 6.23,
                     5.52, 4.90, 4.65, 4.60, 4.73, 4.99, 5.36, 5.91,
                     6.55, 7.45, 8.45, 9.80, 11.30]
        MW_exvebv = np.array(MW_exvebv)
        Rv = 3.08
        MW_axav = MW_exvebv/Rv + 1.0
        ax.plot(1./MW_x, MW_axav, 'o', label='MW Observed')

        ax.set_xscale('log')
        ax.set_yscale('log')

        ax.set_ylim(1e-3,10.)

        ax.set_xlabel(r'$\lambda$ [$\mu$m]')
        ax.set_ylabel(r'$A(x)/A(V)$')

        ax.legend(loc='best')
        plt.show()
    """

    # ['FUV_width', 'NUV_width', 'SIL1_width', 'SIL2_width', 'FIR_width']
    # [0.14696938456699066, 0.051719412985268637, 3.7049250106326643, 8.3409029232536902, 35.355339059327378]

    inputs = ("x",)
    outputs = ("axav",)

    # constant for conversion from Ax/Ab to (more standard) Ax/Av
    AbAv = 1.0 / 3.08 + 1.0

    BKG_amp = Parameter(
        description="BKG term: amplitude", default=165.0 * AbAv, min=0.0
    )
    BKG_lambda = Parameter(description="BKG term: center wavelength", default=0.047)
    BKG_width = Parameter(description="BKG term: width", default=0.452, min=0.0)

    FUV_amp = Parameter(description="FUV term: amplitude", default=14.0 * AbAv, min=0.0)
    FUV_lambda = Parameter(
        description="FUV term: center wavelength", default=0.07, bounds=(0.06, 0.08)
    )
    FUV_b = Parameter(description="FUV term: b coefficient", default=4.0)
    FUV_n = Parameter(description="FUV term: n coefficient", default=6.5)

    NUV_amp = Parameter(
        description="NUV term: amplitude", default=0.045 * AbAv, min=0.0
    )
    NUV_lambda = Parameter(
        description="NUV term: center wavelength", default=0.22, bounds=(0.20, 0.24)
    )
    NUV_width = Parameter(description="NUV term: width", default=0.05, min=0.0)

    SIL1_amp = Parameter(
        description="SIL1 term: amplitude", default=0.002 * AbAv, bounds=(0.0, 0.01)
    )
    SIL1_lambda = Parameter(
        description="SIL1 term: center wavelength", default=9.7, bounds=(7.0, 13.0)
    )
    SIL1_width = Parameter(description="SIL1 term: width", default=2.0, bounds=(0.0, 5.0))

    SIL2_amp = Parameter(
        description="SIL2 term: amplitude", default=0.002 * AbAv, min=0.00
    )
    SIL2_lambda = Parameter(
        description="SIL2 term: center wavelength", default=18.0, bounds=(15.0, 21.0)
    )
    SIL2_width = Parameter(description="SIL2 term: width", default=8.3, min=0.0)

    FIR_amp = Parameter(
        description="FIR term: amplitude", default=0.012 * AbAv, min=0.0
    )
    FIR_lambda = Parameter(
        description="FIR term: center wavelength", default=25.0, bounds=(20.0, 30.0)
    )
    FIR_width = Parameter(description="FIR term: width", default=35.0)

    x_range = [1.0 / 1e3, 1.0 / 1e-3]

    @staticmethod
    def _p92_single_term(in_lambda, amplitude, cen_wave, b, n):
        """
        Function for calculating a single P92 term

        .. math::

           \frac{a}{(\lambda/cen_wave)^n + (cen_wave/\lambda)^n + b}

        when n = 2, this term is equivalent to a Drude profile

        Parameters
        ----------
        in_lambda: vector of floats
           wavelengths in same units as cen_wave

        amplitude: float
           amplitude

        cen_wave: flaot
           central wavelength

        b : float
           b coefficient

        n : float
           n coefficient
        """
        l_norm = in_lambda / cen_wave

        return amplitude / (np.power(l_norm, n) + np.power(l_norm, -1 * n) + b)

    def evaluate(
        self,
        in_x,
        BKG_amp,
        BKG_lambda,
        BKG_width,
        FUV_amp,
        FUV_lambda,
        FUV_b,
        FUV_n,
        NUV_amp,
        NUV_lambda,
        NUV_width,
        SIL1_amp,
        SIL1_lambda,
        SIL1_width,
        SIL2_amp,
        SIL2_lambda,
        SIL2_width,
        FIR_amp,
        FIR_lambda,
        FIR_width,
    ):
        """
        P92 function

        Parameters
        ----------
        in_x: float
           expects either x in units of wavelengths or frequency
           or assumes wavelengths in wavenumbers [1/micron]

           internally wavenumbers are used

        Returns
        -------
        axav: np array (float)
            A(x)/A(V) extinction curve [mag]

        Raises
        ------
        ValueError
           Input x values outside of defined range
        """
        x = _get_x_in_wavenumbers(in_x)

        # check that the wavenumbers are within the defined range
        _test_valid_x_range(x, self.x_range, "P92_mod")

        # compute b from lambda and width
        BKG_b = np.power((BKG_width / BKG_lambda), 2.0) - 2.0
        NUV_b = np.power((NUV_width / NUV_lambda), 2.0) - 2.0
        SIL1_b = np.power((SIL1_width / SIL1_lambda), 2.0) - 2.0
        SIL2_b = np.power((SIL2_width / SIL2_lambda), 2.0) - 2.0
        FIR_b = np.power((FIR_width / FIR_lambda), 2.0) - 2.0

        # calculate the terms
        lam = 1.0 / x
        axav = (
            self._p92_single_term(lam, BKG_amp, BKG_lambda, BKG_b, 2.0)
            + self._p92_single_term(lam, FUV_amp, FUV_lambda, FUV_b, FUV_n)
            + self._p92_single_term(lam, NUV_amp, NUV_lambda, NUV_b, 2.0)
            + self._p92_single_term(lam, SIL1_amp, SIL1_lambda, SIL1_b, 2.0)
            + self._p92_single_term(lam, SIL2_amp, SIL2_lambda, SIL2_b, 2.0)
            + self._p92_single_term(lam, FIR_amp, FIR_lambda, FIR_b, 2.0)
        )

        # return A(x)/A(V)
        return axav

    # use numerical derivaties (need to add analytic)
    fit_deriv = None


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
    # param_dict = dict(zip(model.param_names, model.parameters))
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
    return -0.5 * np.sum(((model(x) - y) / uncs) ** 2)


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


def p92_emcee(
    x,
    y,
    uncs,
    model,
    fit_param_names=None,
    threads=1,
    return_sampler=False,
    nburn=100,
    nsteps=500,
):
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
    nwalkers = 10 * ndim

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

    # print(fit_param_names)
    # print(p0)

    # setting up the walkers to start "near" the inital guess
    p = [p0 * (1 + 1e-4 * np.random.normal(0, 1.0, ndim)) for k in range(nwalkers)]

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
    sampler = emcee.EnsembleSampler(
        nwalkers,
        ndim,
        lnprob,
        threads=threads,
        args=(x, y, uncs, model_copy, fit_param_names),
    )

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
    per_params = [
        (v[1], v[2] - v[1], v[1] - v[0])
        for v in zip(*np.percentile(samples, [16, 50, 84], axis=0))
    ]

    # set the returned model parameters to the best fit values
    for k, val in enumerate(per_params):
        exec("model_copy.{}.value = {}".format(fit_param_names[k], best_params[k]))
        print(fit_param_names[k], best_params[k], val)

    clean_pnames = [pname[:-2] for pname in fit_param_names]
    model_copy.p92_emcee_param_names = clean_pnames
    model_copy.p92_emcee_per_params = per_params

    if return_sampler:
        return (model_copy, sampler)
    else:
        return model_copy


def plot_emcee_results(sampler, fit_param_names, filebase=""):
    """
    Plot the standard triangle and diagnostic walker plots
    """

    # plot the walker chains for all parameters
    nwalkers, nsteps, ndim = sampler.chain.shape
    fig, ax = plt.subplots(ndim, sharex=True, figsize=(13, 13))
    walk_val = np.arange(nsteps)
    for i in range(ndim):
        for k in range(nwalkers):
            ax[i].plot(walk_val, sampler.chain[k, :, i], "-")
            ax[i].set_ylabel(fit_param_names[i])
    fig.savefig("%s_walker_param_values.png" % filebase)
    plt.close(fig)

    # plot the 1D and 2D likelihood functions in a traditional triangle plot
    samples = sampler.chain.reshape((-1, ndim))
    fig = corner.corner(
        samples,
        labels=fit_param_names,
        show_titles=True,
        title_fmt=".3f",
        use_math_text=True,
    )
    fig.savefig("%s_param_triangle.png" % filebase)
    plt.close(fig)


if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("redstarname", help="name of reddened star")
    parser.add_argument("compstarname", help="name of comparision star")
    parser.add_argument(
        "--path",
        help="base path to observed data",
        default="/home/kgordon/Python_git/extstar_data/",
    )
    parser.add_argument("--emcee", help="run EMCEE fit", action="store_true")
    parser.add_argument(
        "--nburn", type=int, default=100, help="# of burn steps in MCMC chain"
    )
    parser.add_argument(
        "--nsteps", type=int, default=500, help="# of steps in MCMC chain"
    )
    parser.add_argument(
        "--threads", type=int, default=1, help="number of threads for EMCEE run"
    )
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    # read in the observed data for both stars
    redstarobs = StarData("DAT_files/%s.dat" % args.redstarname, path=args.path)
    compstarobs = StarData("DAT_files/%s.dat" % args.compstarname, path=args.path)

    # output filebase
    filebase = "fits/%s_%s" % (args.redstarname, args.compstarname)

    # calculate the extinction curve
    extdata = ExtData()
    extdata.calc_elx(redstarobs, compstarobs)

    # get an observed extinction curve to fit
    (wave, y, y_unc) = extdata.get_fitdata(
        ["BAND", "IUE", "IRS"], remove_uvwind_region=True, remove_lya_region=True
    )
    # remove data affected by Ly-alpha absorption/emission
    # gindxs = wave > 1400 * u.Angstrom
    gindxs = wave > 500 * u.Angstrom
    wave = wave[gindxs]
    y = y[gindxs]
    y_unc = y_unc[gindxs]
    # remove units as fitting routines often cannot take numbers with units
    x = wave.to(1.0 / u.micron, equivalencies=u.spectral()).value

    # determine the initial guess at the A(V) values
    #  just use the average at wavelengths > 5
    #  limit as lambda -> inf, E(lamda-V) -> -A(V)
    indxs, = np.where(1.0 / x > 5.0)
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
    p92_init.FUV_lambda_0.fixed = True
    p92_init.FUV_b_0.fixed = True
    p92_init.FUV_n_0.fixed = True
    p92_init.SIL2_lambda_0.fixed = True
    p92_init.SIL2_width_0.fixed = True
    p92_init.FIR_lambda_0.fixed = True
    p92_init.FIR_width_0.fixed = True

    p92_init.Av_1.bounds = [0.1, None]

    # extra for HD 147933 (4 Nov 2019)
    #p92_init.SIL1_amp_0.bounds = [0.005, None]
    #p92_init.BKG_lambda_0.fixed = False
    #p92_init.BKG_width_0.fixed = False
    #p92_init.FUV_lambda_0.fixed = False

    # p92_init.SIL2_amp_0.tied = tie_amps_SIL2_to_SIL1
    # p92_init.FIR_amp_0.tied = tie_amps_FIR_to_SIL1

    # pick the fitter
    fit = LevMarLSQFitter()

    # fit the data to the P92 model using the fitter
    p92_fit = fit(p92_init, x, y, weights=1.0 / y_unc, maxiter=10000, epsilon=0.001)

    clean_pnames = [pname[:-2] for pname in p92_init.param_names]

    p92_fit_params = p92_fit.parameters

    for k, cur_pname in enumerate(clean_pnames):
        print("{:12} {:6.4e}".format(cur_pname, p92_fit_params[k]))

    best_fit_Av = p92_fit.Av_1.value

    p92_best_params = (clean_pnames, p92_fit_params)
    p92_per_params = None
    if args.emcee:
        # run the emcee fitter to get proper fit parameter uncertainties
        fit_param_names = [
            "Av_1",
            "BKG_amp_0",
            "FUV_amp_0",
            "NUV_amp_0",
            "NUV_lambda_0",
            "NUV_width_0",
            "SIL1_amp_0",
            "SIL1_lambda_0",
            "SIL1_width_0",
            "SIL2_amp_0",
            "FIR_amp_0",
        ]
        emcee_results = p92_emcee(
            x,
            y,
            y_unc,
            p92_fit,
            nburn=args.nburn,
            nsteps=args.nsteps,
            threads=args.threads,
            fit_param_names=fit_param_names,
            return_sampler=True,
        )

        p92_fit_emcee, sampler = emcee_results

        clean_pnames_emcee = [pname[:-2] for pname in fit_param_names]

        plot_emcee_results(sampler, clean_pnames_emcee, filebase=filebase)

        p92_per_params = (
            p92_fit_emcee.p92_emcee_param_names,
            p92_fit_emcee.p92_emcee_per_params,
        )
        p92_best_params = (clean_pnames, p92_fit_emcee.parameters)

    # save the extinction curve and fit
    warnings.simplefilter("ignore", category=AstropyWarning)
    out_fname = "fits/%s_%s_ext.fits" % (args.redstarname, args.compstarname)
    extdata.save(
        out_fname, p92_best_params=p92_best_params, p92_per_params=p92_per_params
    )

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
    # fig, ax = plt.subplots(figsize=(12, 8))
    fig, ax = plt.subplots(figsize=(6, 4))

    # subplot
    ax2 = plt.axes([0.60, 0.35, 0.35, 0.35])

    # plot the bands and all spectra for this star
    extdata.plot(ax, color="k", alpha=0.5)
    extdata.plot(ax2, color="k", alpha=0.5)

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
                    exec(
                        "p92_fit_emcee.{}.value = {}".format(
                            pname, sampler.chain[k, j, i]
                        )
                    )
                ax.plot(1.0 / x, p92_fit_emcee(x), "b-", alpha=0.01)
                ax2.plot(1.0 / x, p92_fit_emcee(x), "b-", alpha=0.01)

        # p50 results
        ax.plot(1.0 / x, p92_fit_emcee(x), "b-", label="EMCEE Fits")
        ax2.plot(1.0 / x, p92_fit_emcee(x), "b-")

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
    indxs, = np.where((x > 1.0 / sp_xlim[1]) & (x < 1.0 / sp_xlim[0]))
    ax2.set_ylim(
        min([min(p92_fit(x)[indxs]), -best_fit_Av]) - 0.1, max(p92_fit(x)[indxs]) + 0.1
    )

    # use the whitespace better
    warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")
    fig.tight_layout()

    # plot or save to a file
    outname = "%s_ext" % filebase
    if args.png:
        fig.savefig(outname + ".png")
    elif args.pdf:
        fig.savefig(outname + ".pdf")
    else:
        plt.show()
