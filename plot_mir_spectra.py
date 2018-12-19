#!/usr/bin/env python
#
# Program to plot a list of spectra used for extinction curves
#
# written Dec 2014/Jan 2015 by Karl Gordon (kgordon@stsci.edu)
# based strongly on IDL programs created over the previous 10 years
#
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import argparse

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib

from measure_extinction.stardata import StarData


def plot_mir_set(ax, starnames, extra_off_val=0.0,
                 plam4=True,
                 norm_wave_range=[6., 10.],
                 col_vals=['b', 'g', 'r', 'm', 'c', 'y'],
                 ann_xvals=[35.0, 42.0],
                 fontsize=12):
    """
    Plot a set of spectra
    """

    spec_name = 'IRS'
    for i in range(len(starnames)):
        stardata = StarData('DAT_files/'+starnames[i]+'.dat',
                            path='/home/kgordon/Python_git/extstar_data/')

        if plam4:
            ymult = np.power(stardata.data[spec_name].waves, 4.0)
        else:
            ymult = np.full((len(stardata.data[spec_name].waves)), 1.0)

        # get the value to use for normalization and offset
        norm_indxs = np.where((stardata.data[spec_name].waves
                               >= norm_wave_range[0]) &
                              (stardata.data[spec_name].waves
                               <= norm_wave_range[1]))
        norm_val = 1.0/np.average(
            stardata.data[spec_name].fluxes[norm_indxs]
            * ymult[norm_indxs])
        off_val = extra_off_val + 0.5*i

        # plot the spectroscopic data
        gindxs = np.where(stardata.data[spec_name].npts > 0)
        # max_gwave = max(stardata.data[spec_name].waves[gindxs])
        ax.plot(stardata.data[spec_name].waves[gindxs],
                (stardata.data[spec_name].fluxes[gindxs]*ymult[gindxs]
                 * norm_val + off_val),
                col_vals[i % 6] + '-')

        # annotate the spectra
        # ann_wave_range = np.array([max_gwave-5.0, max_gwave-1.0])
        ann_wave_range = [9.0, 15.0]
        ann_indxs = np.where((stardata.data[spec_name].waves
                              >= ann_wave_range[0]) &
                             (stardata.data[spec_name].waves
                              <= ann_wave_range[1]) &
                             (stardata.data[spec_name].npts > 0))
        ann_val = np.median(stardata.data[spec_name].fluxes[ann_indxs]
                            * ymult[ann_indxs])
        ann_val *= norm_val
        ann_val += off_val + 0.2
        ax.annotate(starnames[i]+' '+stardata.sptype, xy=(ann_xvals[0],
                                                          ann_val),
                    xytext=(ann_xvals[1], ann_val),
                    verticalalignment="center",
                    arrowprops=dict(facecolor=col_vals[i % 6], shrink=0.1),
                    fontsize=0.85*fontsize, rotation=-0.)

        # plot the band fluxes
        if plam4:
            ymult = np.power(stardata.data['BAND'].waves, 4.0)
        else:
            ymult = np.full((len(stardata.data['BAND'].waves)),
                            1.0)
        ax.plot(stardata.data['BAND'].waves,
                stardata.data['BAND'].fluxes*ymult
                * norm_val + off_val, col_vals[i % 6] + 'o')


def ann_set(ax, fontsize,
            bracket_x, bracket_y,
            text_x,
            texta, textb):
    """
    Annotate a set of spectra with text
    """
    ax.plot([bracket_x[0], bracket_x[1], bracket_x[1], bracket_x[0]],
            [bracket_y[0], bracket_y[0], bracket_y[1], bracket_y[1]],
            'k-',
            linewidth=3.)
    text_y = 0.5*(bracket_y[0] + bracket_y[1])
    ax.text(text_x[0], text_y, texta, rotation=270.,
            fontsize=1.2*fontsize,
            horizontalalignment='center', verticalalignment="center")
    ax.text(text_x[1], text_y, textb, rotation=270.,
            fontsize=.9*fontsize,
            horizontalalignment='center', verticalalignment="center")


if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("filelist", help="file with list of stars to plot")
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
    starnames = []
    stardata = []
    for line in file_lines:
        if (line.find('#') != 0) & (len(line) > 0):
            name = line.rstrip()
            starnames.append(name)

    fontsize = 12

    font = {'size': fontsize}

    matplotlib.rc('font', **font)

    matplotlib.rc('lines', linewidth=1)
    matplotlib.rc('axes', linewidth=2)
    matplotlib.rc('xtick.major', width=2)
    matplotlib.rc('xtick.minor', width=2)
    matplotlib.rc('ytick.major', width=2)
    matplotlib.rc('ytick.minor', width=2)

    # fig, ax = pyplot.subplots(nrows=1, ncols=1, figsize=(10, 13))
    fig, ax = pyplot.subplots(nrows=1, ncols=1, figsize=(10, 10))

    # kxrange = [3.0, 130.]
    kxrange = [3.0, 65.]
    # ann_xvals = [41.0, 50.0]
    ann_xvals = [35.0, 42.0]
    spec_name = 'IRS'
    norm_wave_range = [6., 10.]
    # ann_wave_range = [15.0, 18.0]
    col_vals = ['b', 'g', 'r', 'm', 'c', 'y']
    # col_vals = ['k', 'k', 'k', 'k', 'k', 'k']

    plot_mir_set(ax, starnames)

    ax.set_yscale('linear')
    ax.set_ylim(0.5, len(starnames) * 0.5 + 4.0)
    ax.set_xscale('log')
    ax.set_xlim(kxrange)
    ax.set_xlabel('$\lambda$ [$\mu m$]', fontsize=1.3 * fontsize)
    ax.set_ylabel('$\lambda^4 F(\lambda)/F(8 \mu m)$ + offset',
                  fontsize=1.3*fontsize)

    ax.tick_params('both', length=10, width=2, which='major')
    ax.tick_params('both', length=5, width=1, which='minor')

    fig.tight_layout()

    save_str = '_mspec'
    if args.png:
        fig.savefig(args.filelist.replace('.dat', save_str+'.png'))
    elif args.eps:
        fig.savefig(args.filelist.replace('.dat', save_str+'.eps'))
    elif args.pdf:
        fig.savefig(args.filelist.replace('.dat', save_str+'.pdf'))
    else:
        pyplot.show()
