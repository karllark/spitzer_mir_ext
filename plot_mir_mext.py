#!/usr/bin/env python
#
# Program to plot a list of extinction curves
#
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import argparse

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib

from calc_ext import P92_Elv
from dust_extinction.dust_extinction import P92
from measure_extinction.extdata import ExtData

if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("filelist", help="file with list of curves to plot")
    parser.add_argument("--alav", help="plot A(lambda)/A(V)",
                        action="store_true")
    parser.add_argument("--modonly", help="only plot the models",
                        action="store_true")
    parser.add_argument("-p", "--png", help="save figure as a png file",
                        action="store_true")
    parser.add_argument("-e", "--eps", help="save figure as an eps file",
                        action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file",
                        action="store_true")
    args = parser.parse_args()

    filename = args.filelist

    f = open(filename, 'r')
    file_lines = list(f)
    extnames = []
    extdatas = []
    avs = []
    for line in file_lines:
        if (line.find('#') != 0) & (len(line) > 0):
            name = line.rstrip()
            extnames.append(name)
            text = ExtData(filename=name)
            extdatas.append(text)
            avs.append(text.columns['AV'][0])

    fontsize = 18

    font = {'size': fontsize}

    matplotlib.rc('font', **font)

    matplotlib.rc('lines', linewidth=1)
    matplotlib.rc('axes', linewidth=2)
    matplotlib.rc('xtick.major', width=2)
    matplotlib.rc('xtick.minor', width=2)
    matplotlib.rc('ytick.major', width=2)
    matplotlib.rc('ytick.minor', width=2)

    fig, ax = pyplot.subplots(nrows=1, ncols=1, figsize=(10, 13))

    sindxs = np.argsort(avs)

    kxrange = [1.0, 130.]
    ann_xvals = [41.0, 50.0]
    spec_name = 'IRS'
    norm_wave_range = [6., 10.]
    ann_wave_range = [15.0, 18.0]
    col_vals = ['b', 'g', 'r', 'm', 'c', 'y']

    mod_x = 1.0/np.arange(1.0, 40.0, 0.1)
    for i in range(len(extnames)):
        k = sindxs[i]

        # get the value to use for normalization and offset
        norm_indxs = np.where((extdatas[k].waves[spec_name]
                               >= norm_wave_range[0]) &
                              (extdatas[k].waves[spec_name]
                               <= norm_wave_range[1]))

        # plot the extinction curves
        if not args.modonly:
            extdatas[k].plot_ext(ax, color=col_vals[i % 6], alav=args.alav)

        # plot the best fit P92 model
        if args.alav:
            P92_best = P92(BKG_amp=extdatas[k].p92_best_fit['BKG_AMP'],
                           BKG_lambda=extdatas[k].p92_best_fit['BKG_LAMBDA'],
                           BKG_n=extdatas[k].p92_best_fit['BKG_N'],
                           BKG_b=extdatas[k].p92_best_fit['BKG_B'],
                           FUV_amp=extdatas[k].p92_best_fit['FUV_AMP'],
                           FUV_lambda=extdatas[k].p92_best_fit['FUV_LAMBDA'],
                           FUV_n=extdatas[k].p92_best_fit['FUV_N'],
                           FUV_b=extdatas[k].p92_best_fit['FUV_B'],
                           NUV_amp=extdatas[k].p92_best_fit['NUV_AMP'],
                           NUV_lambda=extdatas[k].p92_best_fit['NUV_LAMBDA'],
                           NUV_n=extdatas[k].p92_best_fit['NUV_N'],
                           NUV_b=extdatas[k].p92_best_fit['NUV_B'],
                           SIL1_amp=extdatas[k].p92_best_fit['SIL1_AMP'],
                           SIL1_lambda=extdatas[k].p92_best_fit['SIL1_LAMBDA'],
                           SIL1_n=extdatas[k].p92_best_fit['SIL1_N'],
                           SIL1_b=extdatas[k].p92_best_fit['SIL1_B'],
                           FIR_amp=extdatas[k].p92_best_fit['FIR_AMP'],
                           FIR_lambda=extdatas[k].p92_best_fit['FIR_LAMBDA'],
                           FIR_n=extdatas[k].p92_best_fit['FIR_N'],
                           FIR_b=extdatas[k].p92_best_fit['FIR_B'])
        else:
            P92_best = P92_Elv(BKG_amp_0=extdatas[k].p92_best_fit['BKG_AMP'],
                           BKG_lambda_0=extdatas[k].p92_best_fit['BKG_LAMBDA'],
                           BKG_n_0=extdatas[k].p92_best_fit['BKG_N'],
                           BKG_b_0=extdatas[k].p92_best_fit['BKG_B'],
                           FUV_amp_0=extdatas[k].p92_best_fit['FUV_AMP'],
                           FUV_lambda_0=extdatas[k].p92_best_fit['FUV_LAMBDA'],
                           FUV_n_0=extdatas[k].p92_best_fit['FUV_N'],
                           FUV_b_0=extdatas[k].p92_best_fit['FUV_B'],
                           NUV_amp_0=extdatas[k].p92_best_fit['NUV_AMP'],
                           NUV_lambda_0=extdatas[k].p92_best_fit['NUV_LAMBDA'],
                           NUV_n_0=extdatas[k].p92_best_fit['NUV_N'],
                           NUV_b_0=extdatas[k].p92_best_fit['NUV_B'],
                           SIL1_amp_0=extdatas[k].p92_best_fit['SIL1_AMP'],
                           SIL1_lambda_0=extdatas[k].p92_best_fit['SIL1_LAMBDA'],
                           SIL1_n_0=extdatas[k].p92_best_fit['SIL1_N'],
                           SIL1_b_0=extdatas[k].p92_best_fit['SIL1_B'],
                           FIR_amp_0=extdatas[k].p92_best_fit['FIR_AMP'],
                           FIR_lambda_0=extdatas[k].p92_best_fit['FIR_LAMBDA'],
                           FIR_n_0=extdatas[k].p92_best_fit['FIR_N'],
                           FIR_b_0=extdatas[k].p92_best_fit['FIR_B'],
                           Av_1=extdatas[k].columns['AV'][0])

        ax.plot(1.0/mod_x, P92_best(mod_x), '--',
                color=col_vals[i % 6], alpha=0.5)

        # annotate the spectra
        #ann_wave_range = np.array([max_gwave-5.0, max_gwave-1.0])
        #ann_indxs = np.where((extdata[k].data[spec_name].waves
        #                      >= ann_wave_range[0]) &
        #                     (extdata[k].data[spec_name].waves
        #                      <= ann_wave_range[1]) &
        #                     (extdata[k].data[spec_name].npts > 0))
        #ann_val = np.median(extdata[k].data[spec_name].fluxes[ann_indxs]
        #                    *ymult[ann_indxs])
        #ann_val *= norm_val
        #ann_val += off_val
        #ax.annotate(starnames[k]+'/'+extdata[k].sptype, xy=(ann_xvals[0],
        #                                                     ann_val),
        #            xytext=(ann_xvals[1], ann_val),
        #            verticalalignment="center",
        #            arrowprops=dict(facecolor=col_vals[i%6], shrink=0.1),
        #            fontsize=0.85*fontsize, rotation=-0.)

    ax.set_yscale('linear')
    ax.set_xscale('log')
    ax.set_xlim(kxrange)
    if args.alav:
        ax.set_ylim(-0.05, 0.2)
        ax.set_ylabel('$A(\lambda)/A(V)$',
                      fontsize=1.3*fontsize)
        # Milky Way observed extinction from
        # Rieke & Lebofsky (1985)
        MW_x = 1.0/np.array([13.0, 12.5, 12.0, 11.5, 11.0, 10.5,
                            10.0,  9.5,  9.0,  8.5,  8.0,  4.8,
                            3.5, 2.22, 1.65, 1.25])
        MW_axav = np.array([0.027, 0.030, 0.037, 0.047, 0.060, 0.074,
                            0.083, 0.087, 0.074, 0.043, 0.020, 0.023,
                            0.058, 0.112, 0.175, 0.282])
        ax.plot(1.0/MW_x, MW_axav, 'ko')

    else:
        ax.set_ylim(-6.0, -1.0)
        ax.set_ylabel('$E(\lambda - V)$',
                      fontsize=1.3*fontsize)

    ax.tick_params('both', length=10, width=2, which='major')
    ax.tick_params('both', length=5, width=1, which='minor')

    fig.tight_layout()

    save_str = '_mext'
    if args.png:
        fig.savefig(args.filelist.replace('.dat', save_str+'.png'))
    elif args.eps:
        fig.savefig(args.filelist.replace('.dat', save_str+'.eps'))
    elif args.pdf:
        fig.savefig(args.filelist.replace('.dat', save_str+'.pdf'))
    else:
        pyplot.show()
