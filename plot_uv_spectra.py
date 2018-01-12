#!/usr/bin/env python2.7
#
# Program to plot a list of spectra used for extinction curves
#
# written Dec 2014/Jan 2015 by Karl Gordon (kgordon@stsci.edu)
# based strongly on IDL programs created over the previous 10 years
#
from __future__ import print_function
import argparse

import numpy as np
from astropy.io import fits as pyfits
import matplotlib.pyplot as pyplot
import matplotlib

from getstardata import *

if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("filelist", help="file with list of stars to plot")
    parser.add_argument("-p", "--png", help="save figure as a png file",
                        action="store_true")
    parser.add_argument("-e", "--eps", help="save figure as an eps file",
                        action="store_true")
    args = parser.parse_args()

    filename = args.filelist

    f = open(filename, 'r')
    file_lines = list(f)
    starnames = []
    stardata = []
    bvcol = []
    for line in file_lines:
        if (line.find('#') != 0) & (len(line) > 0):
            name = line.rstrip()
            starnames.append(name)
            tstar = StarData('DAT_files/'+name+'.dat', path='/home/kgordon/Dust/Ext/')
            stardata.append(tstar)

    fontsize = 18

    font = {'size'   : fontsize}

    matplotlib.rc('font', **font)

    matplotlib.rc('lines', linewidth=1)
    matplotlib.rc('axes', linewidth=2)
    matplotlib.rc('xtick.major', width=2)
    matplotlib.rc('ytick.major', width=2)
            
    fig, ax = pyplot.subplots(nrows=1,ncols=1, figsize=(10,12))

    #sindxs = np.argsort(bvcol)

    norm_wave_range = 1.0/np.array([7.5,7.0])
    ann_wave_range = 1.0/np.array([7.65,7.45])
    col_vals = ['b','g','r','m','c','y']
    for i in range(len(starnames)):
        #k = sindxs[len(starnames) - i - 1]
        k = i
        #norm_indxs = np.where((stardata[k].data['BANDS'].flat_bands_waves >= norm_wave_range[0]) &
        #                      (stardata[k].data['BANDS'].flat_bands_waves <= norm_wave_range[1]))
        #norm_val = 1.0/np.average(stardata[k].data['BANDS'].flat_bands_fluxes[norm_indxs])

        if 'STIS' in stardata[k].data.keys():

            norm_indxs = np.where((stardata[k].data['STIS'].waves >= norm_wave_range[0]) &
                                  (stardata[k].data['STIS'].waves <= norm_wave_range[1]))
            norm_val = 1.0/np.average(stardata[k].data['STIS'].flux[norm_indxs])
            
            norm_val *= 6.0**i
            
            ax.plot(1.0/stardata[k].data['BANDS'].flat_bands_waves,stardata[k].data['BANDS'].flat_bands_fluxes*norm_val,col_vals[i%6]+'o')

            gindxs = np.where(stardata[k].data['STIS'].npts > 0)
            ax.plot(1.0/stardata[k].data['STIS'].waves[gindxs],stardata[k].data['STIS'].flux[gindxs]*norm_val,col_vals[i%6]+'-')

            
            ann_indxs = np.where((stardata[k].data['STIS'].waves >= ann_wave_range[0]) &
                                 (stardata[k].data['STIS'].waves <= ann_wave_range[1]) &
                                 (stardata[k].data['STIS'].npts > 0))
            ann_val = np.median(stardata[k].data['STIS'].flux[ann_indxs])
            ann_val *= norm_val
            ax.annotate(starnames[k], xy=(8.8, ann_val), xytext=(9.3, ann_val), 
                        verticalalignment="center",
                        arrowprops=dict(facecolor=col_vals[i%6], shrink=0.1), 
                        fontsize=0.85*fontsize, rotation=-45.)

    #ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(1e-2,7e17)
    #ax.set_ylim(1e-2,7e14)
    ax.set_xlim(3.0,10.35)
    ax.set_xlabel('$1/\lambda$ [$\mu m^{-1}$]',fontsize=1.3*fontsize)
    ax.set_ylabel('$F(\lambda)/F(3.5)$ + offset',fontsize=1.3*fontsize)
    
    fig.tight_layout()

    save_str = '_mspec'
    if args.png:
        fig.savefig(string.replace(args.filelist,'.dat',save_str+'.png'))
    elif args.eps:
        fig.savefig(string.replace(args.filelist,'.dat',save_str+'.eps'))
    else:
        pyplot.show()

