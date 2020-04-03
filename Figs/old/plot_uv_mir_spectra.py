#!/usr/bin/env python

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import argparse

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

from measure_extinction.stardata import StarData

if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("filelist", help="file with list of stars to plot")
    parser.add_argument("--png", help="save figure as a png file",
                        action="store_true")
    parser.add_argument("--eps", help="save figure as an eps file",
                        action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file",
                        action="store_true")
    args = parser.parse_args()

    filename = args.filelist

    f = open(filename, 'r')
    file_lines = list(f)
    starnames = []
    stardata = []
    for line in file_lines:
        if (line.find('#') != 0) & (len(line) > 0):
            name = line.rstrip()
            starnames.append(name)
            tstar = StarData('DAT_files/'+name+'.dat',
                             path='/home/kgordon/Dust/Ext/')
            stardata.append(tstar)

    fontsize = 14

    font = {'size': fontsize}
    mpl.rc('font', **font)
    mpl.rc('lines', linewidth=1)
    mpl.rc('axes', linewidth=2)
    mpl.rc('xtick.major', width=2)
    mpl.rc('xtick.minor', width=2)
    mpl.rc('ytick.major', width=2)
    mpl.rc('ytick.minor', width=2)

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(16, 6))

    uvoffval = 1.0
    miroffval = 0.5
    pcolors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    for k in range(len(starnames)):

        # plot UV data
        stardata[k].plot_obs(ax[0], pcolor=pcolors[k % len(pcolors)],
                             yoffset=k*uvoffval,
                             norm_wave_range=[0.28, 0.3])

        # plot MIR data
        plotvals = stardata[k].plot_obs(ax[1],
                                        pcolor=pcolors[k % len(pcolors)],
                                        mlam4=True,
                                        yoffset=k*miroffval,
                                        norm_wave_range=[10., 20.],
                                        legend_key='IRS',
                                        fontsize=fontsize)

    # finish the setup of the UV plot
    cax = ax[0]
    cax.set_yscale('log')
    cax.set_ylim(0.01, len(starnames)*uvoffval+10.)
    cax.set_xscale('linear')
    cax.set_xlim(0.1, 0.33)
    cax.set_xlabel('$\lambda$ [$\mu m$]', fontsize=1.3*fontsize)
    cax.set_ylabel('$F(\lambda)/F(\lambda_n)$ + offset',
                   fontsize=1.3*fontsize)
    cax.tick_params('both', length=10, width=2, which='major')
    cax.tick_params('both', length=5, width=1, which='minor')
    # finish the setup of the MIR plot
    cax = ax[1]
    cax.set_yscale('linear')
    cax.set_ylim(0.5, len(starnames)*miroffval+1.5)
    cax.set_xscale('linear')
    cax.set_xlim(3.0, 35.)
    cax.set_xlabel('$\lambda$ [$\mu m$]', fontsize=1.3*fontsize)
    cax.set_ylabel('$\lambda^4 F(\lambda)/F(\lambda_n)$ + offset',
                   fontsize=1.3*fontsize)
    cax.tick_params('both', length=10, width=2, which='major')
    cax.tick_params('both', length=5, width=1, which='minor')
    cax.yaxis.tick_right()
    cax.yaxis.set_label_position("right")

    ax[1].legend()

    fig.tight_layout()

    save_str = '_mspec'
    if args.png:
        fig.savefig(args.filelist.replace('.dat', save_str+'.png'))
    elif args.eps:
        fig.savefig(args.filelist.replace('.dat', save_str+'.eps'))
    elif args.pdf:
        fig.savefig(args.filelist.replace('.dat', save_str+'.pdf'))
    else:
        plt.show()
