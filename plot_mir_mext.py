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
from measure_extinction.extdata import ExtData

if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("filelist", help="file with list of curves to plot")
    parser.add_argument("--alav", help="plot A(lambda)/A(V)",
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

    kxrange = [3.0, 130.]
    ann_xvals = [41.0, 50.0]
    spec_name = 'IRS'
    norm_wave_range = [6., 10.]
    ann_wave_range = [15.0, 18.0]
    col_vals = ['b', 'g', 'r', 'm', 'c', 'y']

    mod_x = np.arange(3.0, 40.0, 0.1)
    for i in range(len(extnames)):
        k = sindxs[i]

        # get the value to use for normalization and offset
        norm_indxs = np.where((extdatas[k].waves[spec_name]
                               >= norm_wave_range[0]) &
                              (extdatas[k].waves[spec_name]
                               <= norm_wave_range[1]))

        # plot the extinction curves
        extdatas[k].plot_ext(ax, color=col_vals[i % 6], alav=args.alav)

        # plot the best fit P92 model
        P92_best = P92_Elv(BKG_amp_0=extdatas[k].p92_best_fit['BKG_AMP'],
                           Av_1=extdatas[k].columns['AV'][0])
        print(P92_best(mod_x))
        ax.plot(1.0/mod_x, P92_best(mod_x), '-', color=col_vals[i % 6])

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
