import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import argparse

import astropy.units as u
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.modeling.fitting import SherpaFitter

from dust_extinction.dust_extinction import P92

from extdata import ExtData

from astropy.modeling import (Fittable1DModel, Parameter)

class AlAvToElv(Fittable1DModel):
    """
    Model to convert from A(lambda)/A(V) to E(lambda-V)

    Paramters
    ---------
    Av : float
      dust column in A(V) [mag]
    """
    inputs = ('alav',)
    outputs = ('elv',)

    Av = Parameter(description="A(V)",
                   default=1.0, min=0.0)

    @staticmethod
    def evaluate(alav, Av):
        """
        AlAvToElv function

        Paramters
        ---------
        alav : np array (float)
           E(lambda-V)/E(B-V) values

        Returns
        -------
        elv : np array (float)
           E(lamda - V)
        """
        return (alav - 1.0)*Av

    # use numerical derivaties (need to add analytic)
    fit_deriv = None
    
class P92_Av(P92 | AlAvToElv):
    """ evalute P92 on E(l-V) data including solving for A(V)"""
    
if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="file with the extinction curve to fit")
    parser.add_argument("-p", "--png", help="save figure as a png file",
                        action="store_true")
    parser.add_argument("-e", "--eps", help="save figure as an eps file",
                        action="store_true")
    parser.add_argument("--path", help="path for the extinction curves")
    args = parser.parse_args()

    if args.path:
        locpath = args.path + '/'
    else:
        locpath = ''
    
    # plot info
    fontsize = 18
    font = {'size'   : fontsize}
    matplotlib.rc('font', **font)
    matplotlib.rc('lines', linewidth=1)
    matplotlib.rc('axes', linewidth=2)
    matplotlib.rc('xtick.major', width=2)
    matplotlib.rc('ytick.major', width=2)

    # read in the observed E(l-V) extinction curve
    obsext = ExtData()
    obsext.read_ext_data(locpath+args.file)

    # get the photometric band data
    xdata = np.array(obsext.ext_waves['BANDS'],dtype=np.float64)
    ydata = np.array(obsext.ext_curve['BANDS'],dtype=np.float64)
    ydata_unc = np.array(obsext.ext_curve_uncs['BANDS'],dtype=np.float64)
    
    # get the stis data and concatenate with existing data
    stis_indxs, = np.where(obsext.ext_curve_uncs['STIS'] > 0.)
    xdata = np.concatenate((xdata,obsext.ext_waves['STIS'][stis_indxs]))
    ydata = np.concatenate((ydata,obsext.ext_curve['STIS'][stis_indxs]))
    ydata_unc = np.concatenate((ydata_unc,
                                obsext.ext_curve_uncs['STIS'][stis_indxs]))

    # get the irs data and concatenate with existing data
    irs_indxs, = np.where(obsext.ext_curve_uncs['IRS'] != 0.)
    xdata = np.concatenate((xdata,obsext.ext_waves['IRS'][irs_indxs]))
    ydata = np.concatenate((ydata,obsext.ext_curve['IRS'][irs_indxs]))
    ydata_unc = np.concatenate((ydata_unc,
                                obsext.ext_curve_uncs['IRS'][irs_indxs]))

    # sort data
    sindxs = np.argsort(xdata)
    xdata = xdata[sindxs]
    ydata = ydata[sindxs]
    ydata_unc = ydata_unc[sindxs]

    # check for zero uncertainty points
    
    
    # add units to the data
    x = xdata * u.micron
    y = ydata
    y_unc = ydata_unc
    
    # determine the initial guess at the A(V) values
    #  just use the average at wavelengths > 5
    #  limit as lambda -> inf, E(lamda-V) -> -A(V)
    indxs, = np.where(xdata > 5.0)
    av_guess = -1.0*np.average(ydata[indxs])

    # initialize the Pei (1992) plus A(V) compositemodel
    #    a few tweaks to the starting parameters helps find the solution
    p92_init = P92_Av(BKG_amp_0=165.,
                      FUV_amp_0=100.0,FUV_lambda_0=0.06,
                      Av_1=av_guess)

    # print the initial model info
    print('P92+A(V) model info')
    print(p92_init)
    
    # toggle fit parameters to be fixed
    #p92_init.FIR_amp_0.fixed = True

    # specify and run the fit
    fit = LevMarLSQFitter()
    p92_fit = fit(p92_init, 1./xdata, y, #maxiter=200, acc=1e-10,
                  weights=1.0/y_unc)
    
    # message from the fitter (maxfev is maxiter)
    print('fit message')
    print(fit.fit_info['message'])
    
    # print the initial and final fit parameter values
    #print(p92_fit.__dict__.keys())
    print('initital parameters')
    print(p92_init._parameters)
    print('final parameters')
    print(p92_fit._parameters)

    # run a different fitter
    #sfit = SherpaFitter(statistic='chi2', optimizer='moncar',
    #                    estmethod='confidence')    
    #p92_sfit = sfit(p92_init, 1./xdata, y)
    
    # setup the plot
    fig, ax = plt.subplots(figsize=(12,8))
    # create the inset subplot
    ax2 = plt.axes([.50, .35, .35, .35])

    # plot the observed data
    ax.plot(obsext.ext_waves['BANDS'],obsext.ext_curve['BANDS'],'bo')
    ax2.plot(obsext.ext_waves['BANDS'],obsext.ext_curve['BANDS'],'bo')

    ax.plot(obsext.ext_waves['STIS'][stis_indxs],
            obsext.ext_curve['STIS'][stis_indxs])
    ax2.plot(obsext.ext_waves['STIS'][stis_indxs],
             obsext.ext_curve['STIS'][stis_indxs])
    ax.plot(obsext.ext_waves['IRS'][irs_indxs],
            obsext.ext_curve['IRS'][irs_indxs])
    ax2.plot(obsext.ext_waves['IRS'][irs_indxs],
             obsext.ext_curve['IRS'][irs_indxs])

    # generate x values to provide smooth model curves 
    lam = np.logspace(-3.0, 3.0, num=1000)
    x = (1.0/lam)/u.micron
    ax.plot(1./x, p92_fit(x), label='fit')
    ax2.plot(1./x, p92_fit(x), label='fit')
    ax.plot(1./x, p92_init(x), label='init')

    # plot paramters to make nice figures
    ax.set_xscale('log')
    ax.set_xlim(0.1,40.0)
    ax.set_ylim(-av_guess*1.1,15.0)

    ax.set_xlabel(r'$\lambda$ [$\mu m$]',fontsize=1.3*fontsize)
    ax.set_ylabel(r'$E(\lambda - V)$',fontsize=1.3*fontsize)

    ax2.set_xlim(3.0,40.0)
    ax2.set_ylim(-p92_fit.Av_1.value-0.1,-p92_fit.Av_1.value+0.5)
    
    ax.legend(loc='best')
    plt.show()
    
