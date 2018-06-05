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

from astropy.table import Table

from calc_ext import P92_Elv
from dust_extinction.shapes import P92
from measure_extinction.extdata import (ExtData, AverageExtData)

if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("filelist", help="file with list of curves to plot")
    parser.add_argument("--rebin_fac", type=int, default=None,
                        help="rebin factor for spectra")
    parser.add_argument("--alav", help="plot A(lambda)/A(V)",
                        action="store_true")
    parser.add_argument("--ave", help="plot the average",
                        action="store_true")
    parser.add_argument("--models", help="plot the best fit models",
                        action="store_true")
    parser.add_argument("--modonly", help="only plot the best fit models",
                        action="store_true")
    parser.add_argument("--prevobs", help="plot previous observations",
                        action="store_true")
    parser.add_argument("--dg_models", help="plot dust grain models",
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

    if args.alav:
        figsize = (10, 6)
    else:
        figsize = (10, 12)
    fig, ax = pyplot.subplots(nrows=1, ncols=1, figsize=figsize)

    sindxs = np.argsort(avs)

    kxrange = [1.0, 40]
    ann_xvals = [41.0, 50.0]
    spec_name = 'IRS'
    norm_wave_range = [6., 10.]
    ann_wave_range = [15.0, 18.0]
    col_vals = ['b', 'g', 'r', 'm', 'c', 'y']
    lin_vals = ['--', ':', '-.']

    mod_x = 1.0/np.arange(1.0, 40.0, 0.1)
    for i in range(len(extnames)):
        k = sindxs[i]

        # plot the extinction curves
        if not args.modonly:
            extdatas[k].plot_ext(ax, color=col_vals[i % 6],
                                 alav=args.alav,
                                 alpha=0.5,
                                 rebin_fac=args.rebin_fac,
                                 fontsize=fontsize)
                                 # legend_key='IRS')
                                 # annotate_key='IRS',

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

        if args.models:
            if args.alav:
                ltext = None
            else:
                ltext = extdatas[k].red_file.replace('DAT_files/','')
                ltext = ltext.replace('.dat','')
            ax.plot(1.0/mod_x, P92_best(mod_x), lin_vals[i%3],
                    color=col_vals[i % 6], alpha=0.5,
                    label=ltext)

    ax.set_yscale('linear')
    ax.set_xscale('log')
    ax.set_xlim(kxrange)
    if args.alav:
        ax.set_ylim(0.0, 0.25)
        ax.set_ylabel('$A(\lambda)/A(V)$',
                      fontsize=1.3*fontsize)
        if args.prevobs:
            # Milky Way observed extinction from
            # Rieke & Lebofsky (1985)
            MW_x = 1.0/np.array([13.0, 12.5, 12.0, 11.5, 11.0, 10.5,
                                10.0,  9.5,  9.0,  8.5,  8.0,  4.8,
                                3.5, 2.22, 1.65, 1.25])
            MW_axav = np.array([0.027, 0.030, 0.037, 0.047, 0.060, 0.074,
                                0.083, 0.087, 0.074, 0.043, 0.020, 0.023,
                                0.058, 0.112, 0.175, 0.282])
            ax.plot(1.0/MW_x, MW_axav, 'bo',
                    label='GalCenter; Rieke & Lebofsky (1985)')

            remy_wave = np.array([1.24,1.664,2.164,3.545,4.442,5.675,7.760])
            remy_alak = np.array([2.50,1.55,1.00,0.56,0.43,0.43,0.43])
            remy_alak_unc = np.array([0.15,0.08,0.0,0.06,0.08,0.10,0.10])
            avak = 8.5
            remy_alav = remy_alak/avak
            remy_alav_unc = remy_alak_unc/avak
            # ax.errorbar(remy_wave, remy_alav, yerr=remy_alav_unc,
            #             fmt='go', label='Indebetouw et al. (2005)')
            ax.plot(remy_wave, remy_alav, 'go',
                    label='GalPlane; Indebetouw et al. (2005)')

            lutz_wave = np.array([2.622, 2.748, 2.852, 3.017, 3.282, 3.742,
                                  3.996, 4.347, 5.097, 5.865, 6.749, 7.411,
                                  8.609, 12.28, 18.72])
            lutz_alav = np.array([0.07519, 0.06481, 0.07513, 0.08142, 0.06922,
                                  0.05500, 0.05070, 0.05205, 0.04859, 0.05054,
                                  0.04532, 0.04349, 0.07962, 0.05445, 0.05499])
            # ax.plot(lutz_wave, lutz_alav, 'mo', label='GalCenter; Lutz (1999)')

            a = Table.read('data/fritz11_galcenter.dat',
                           format='ascii.commented_header')
            ax.plot(a['wave'], 0.12*a['ext']/2.62, 'mo',
                    label='GalCenter; Fritz et al. (2011)')

            a = Table.read('data/pixie_dust_chiar_2005_modified.dat',
                           format='ascii.commented_header')
            # table in A(l)/A(K) units
            #ax.plot(a['wave'], 0.12*a['galcen'], 'co',
            #        label='Chiar & Tielens (2005)')
            ax.plot(a['wave'], 0.12*a['local'], 'co',
                    label='GalCenter; Chiar & Tielens (2005)')

            if args.dg_models:
                a = Table.read('data/kext_albedo_WD_MW_3.1_60_D03.all_modified',
                               format='ascii.commented_header')
                ax.plot(a['lambda'], a['C_ext/H']/4.802e-22, 'k--',
                        label='MW R(V)=3.1 (Draine 2003)', alpha=0.5)

                b = Table.read('data/kext_albedo_WD_MW_5.5A_30_D03.all_modified',
                               format='ascii.commented_header')
                ax.plot(b['lambda'], b['C_ext/H']/6.622E-22, 'k:',
                        label='MW R(V)=5.5; sizedist=A (Draine 2003)',
                        alpha=0.5)

                b = Table.read('data/kext_albedo_WD_MW_5.5B_30.dat_modified',
                               format='ascii.commented_header')
                ax.plot(b['lambda'], b['C_ext/H']/4.789E-22, 'k-.',
                        label='MW R(V)=5.5; sizedist=B (Weingartner & Draine 2001)',
                        alpha=0.5)


        # get the average extinction curve
        if args.ave:
            ave_extdata = AverageExtData(extdatas, alav=True)
            ave_extdata.plot_ext(ax, color='k',
                                 fontsize=fontsize,
                                 legend_key='IRS',
                                 legend_label='Average (this work)')

        ax.legend(fontsize=12)
    else:
        ax.set_xlim(1.0,40.)
        ax.set_ylim(-6, -0.5)
        ax.set_ylabel('$E(\lambda - V)$',
                      fontsize=1.3*fontsize)
        ax.legend(fontsize=12, ncol=4)

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
