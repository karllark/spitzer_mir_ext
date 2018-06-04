#!/usr/bin/env python
#
# Program to plot the silicate strength versus 2175 A bump strength
#
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import argparse

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib

from measure_extinction.extdata import ExtData

if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("filelist", help="file with list of curves to plot")
    parser.add_argument("--png", help="save figure as a png file",
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
            text = ExtData(filename='fits/%s' % name)
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

    figsize = (8, 8)
    fig, ax = pyplot.subplots(nrows=1, ncols=1, figsize=figsize)

    #sil_amp = extdatas[:].p92_best_fit['SIL1_AMP']
    #print(sil_amp)
    #exit()

    n_ext = len(extdatas)
    sil_amp = np.full((n_ext), 0.0)
    nuv_amp = np.full((n_ext), 0.0)
    for k, cext in enumerate(extdatas):
        sil_amp[k] = cext.p92_best_fit['SIL1_AMP']
        nuv_amp[k] = cext.p92_best_fit['NUV_AMP']

    ax.plot(sil_amp, nuv_amp, 'ko')

    ax.set_yscale('linear')
    ax.set_xscale('linear')
    ax.set_xlabel(r'P92 Silicate 10 $\mu$m amplitude')
    ax.set_ylabel(r'P92 2175 $\AA$ bump amplitude')
    # ax.set_xlim(kxrange)

    ax.tick_params('both', length=10, width=2, which='major')
    ax.tick_params('both', length=5, width=1, which='minor')

    fig.tight_layout()

    save_str = '_sil1_nuv'
    if args.png:
        fig.savefig(args.filelist.replace('.dat', save_str+'.png'))
    elif args.pdf:
        fig.savefig(args.filelist.replace('.dat', save_str+'.pdf'))
    else:
        pyplot.show()
