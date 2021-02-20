#!/usr/bin/env python
#
# Program to plot a list of extinction curves
#
from __future__ import absolute_import, division, print_function, unicode_literals
import argparse
from matplotlib.ticker import ScalarFormatter

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib

from astropy.table import Table

from dust_extinction.averages import (
    RL85_MWGC,
    I05_MWAvg,
    CT06_MWGC,
    CT06_MWLoc,
    F11_MWGC,
)

from measure_extinction.extdata import ExtData, AverageExtData

if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("filelist", help="file with list of curves to plot")
    parser.add_argument(
        "--rebin_fac", type=int, default=None, help="rebin factor for spectra"
    )
    parser.add_argument("--alav", help="plot A(lambda)/A(V)", action="store_true")
    parser.add_argument("--rel_band", help="Band to use for normalization",
                        default="V")
    parser.add_argument("--ave", help="plot the average", action="store_true")
    parser.add_argument(
        "--models", help="plot the best fit models", action="store_true"
    )
    parser.add_argument(
        "--modonly", help="only plot the best fit models", action="store_true"
    )
    parser.add_argument(
        "--prevobs", help="plot previous observations", action="store_true"
    )
    parser.add_argument(
        "--dg_models", help="plot dust grain models", action="store_true"
    )
    parser.add_argument(
        "-p", "--png", help="save figure as a png file", action="store_true"
    )
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    filename = args.filelist

    f = open(filename, "r")
    file_lines = list(f)
    extnames = []
    extdatas = []
    avs = []
    for line in file_lines:
        if (line.find("#") != 0) & (len(line) > 0):
            name = line.rstrip()
            extnames.append(name)
            text = ExtData(filename="fits/%s" % name)
            extdatas.append(text)
            avs.append(text.columns["AV"][0])

    fontsize = 18

    font = {"size": fontsize}

    matplotlib.rc("font", **font)

    matplotlib.rc("lines", linewidth=1)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("xtick.minor", width=2)
    matplotlib.rc("ytick.major", width=2)
    matplotlib.rc("ytick.minor", width=2)

    if args.alav:
        figsize = (10, 8)
    else:
        figsize = (10, 10)
    fig, ax = pyplot.subplots(nrows=1, ncols=1, figsize=figsize)

    sindxs = np.argsort(avs)

    kxrange = [1.0, 100.0]
    ann_xvals = [41.0, 50.0]
    spec_name = "IRS"
    norm_wave_range = [6.0, 10.0]
    ann_wave_range = [15.0, 18.0]
    col_vals = ["b", "g", "c"]  # , "r", "m", "c", "y"]
    lin_vals = ["--", ":", "-."]
    n_cols = len(col_vals)

    mod_x = 1.0 / np.arange(1.0, 40.0, 0.1)
    for i in range(len(extnames)):
        k = sindxs[i]
        # k = i

        # plot the extinction curves
        if not args.modonly:
            extdatas[k].plot(
                ax,
                color=col_vals[i % n_cols],
                alax=args.alav,
                alpha=0.5,
                rebin_fac=args.rebin_fac,
                fontsize=fontsize,
                plot_rel_band=args.rel_band,
            )

            # label the curves
            ltext = extdatas[k].red_file.replace("DAT_files/", "")
            ltext = ltext.replace(".dat", "")

            (wave, y, y_unc) = extdatas[k].get_fitdata(["BAND", "IRS"])
            (indxs,) = np.where(wave.value > 5.0)
            av_guess = -1.0 * np.average(y[indxs])

            if not args.alav:
                if ltext == "vicyg8a":
                    av_guess += 0.1
                elif ltext == "hd112272":
                    av_guess += 0.12
                elif ltext == "hd281159":
                    av_guess += 0.22
                elif ltext == "hd147701":
                    av_guess += 0.17
                elif ltext == "hd229238":
                    av_guess += 0.1
                elif ltext == "hd029309":
                    av_guess += 0.025
                ax.text(
                    40.0,
                    -1.0 * av_guess,
                    ltext,
                    color=col_vals[i % n_cols],
                    fontsize=0.8 * fontsize,
                )

        if args.models:
            # plot the best fit P92 model
            if args.alav:
                P92_best = P92(
                    BKG_amp=extdatas[k].p92_best_fit["BKG_AMP"],
                    BKG_lambda=extdatas[k].p92_best_fit["BKG_LAMBDA"],
                    BKG_width=extdatas[k].p92_best_fit["BKG_WIDTH"],
                    FUV_amp=extdatas[k].p92_best_fit["FUV_AMP"],
                    FUV_lambda=extdatas[k].p92_best_fit["FUV_LAMBDA"],
                    FUV_n=extdatas[k].p92_best_fit["FUV_N"],
                    FUV_b=extdatas[k].p92_best_fit["FUV_B"],
                    NUV_amp=extdatas[k].p92_best_fit["NUV_AMP"],
                    NUV_lambda=extdatas[k].p92_best_fit["NUV_LAMBDA"],
                    NUV_width=extdatas[k].p92_best_fit["NUV_WIDTH"],
                    SIL1_amp=extdatas[k].p92_best_fit["SIL1_AMP"],
                    SIL1_lambda=extdatas[k].p92_best_fit["SIL1_LAMBDA"],
                    SIL1_width=extdatas[k].p92_best_fit["SIL1_WIDTH"],
                    SIL2_amp=extdatas[k].p92_best_fit["SIL2_AMP"],
                    SIL2_lambda=extdatas[k].p92_best_fit["SIL2_LAMBDA"],
                    SIL2_width=extdatas[k].p92_best_fit["SIL2_WIDTH"],
                    FIR_amp=extdatas[k].p92_best_fit["FIR_AMP"],
                    FIR_lambda=extdatas[k].p92_best_fit["FIR_LAMBDA"],
                    FIR_width=extdatas[k].p92_best_fit["FIR_WIDTH"],
                )
            else:
                P92_best = None

            if args.alav:
                ltext = None
            else:
                ltext = extdatas[k].red_file.replace("DAT_files/", "")
                ltext = ltext.replace(".dat", "")
            ax.plot(
                1.0 / mod_x,
                P92_best(mod_x),
                lin_vals[i % 3],
                color=col_vals[i % n_cols],
                alpha=0.5,
                label=ltext,
            )

    ax.set_yscale("linear")
    ax.set_xscale("log")
    # ax.set_xlim(kxrange)
    if args.alav:
        # ax.set_ylim(-0.05, 0.4)
        # ax.set_xlim(1.0, 40.0)
        ax.set_ylabel(rf"$A(\lambda)/A({args.rel_band})$", fontsize=1.3 * fontsize)
        if args.prevobs:
            litmods = [RL85_MWGC(), I05_MWAvg(), CT06_MWGC(), CT06_MWLoc(), F11_MWGC()]
            litdesc = [
                "GalCenter: Rieke & Lebofsky (1985)",
                "GalPlane: Indebetouw et al. (2005)",
                "GalCenter: Chiar & Tielens (2006)",
                "Local: Chiar & Tielens (2006)",
                "GalCenter: Fritz et al. (2011)",
            ]
            litfmt = ["bs", "gP", "c--", "c:", "m^"]
            for k, cmod in enumerate(litmods):
                lit_wave = 1.0 / cmod.obsdata_x
                lit_axav = cmod.obsdata_axav

                ax.plot(
                    lit_wave, lit_axav, litfmt[k], alpha=0.25, label=litdesc[k], lw=3,
                )

        if args.dg_models:
            a = Table.read(
                "data/old/kext_albedo_WD_MW_3.1_60_D03.all_modified",
                format="ascii.commented_header",
            )
            ax.plot(
                a["lambda"],
                a["C_ext/H"] / 4.802e-22,
                "k--",
                label="MW R(V)=3.1 (Draine 2003)",
                alpha=0.5,
            )

            b = Table.read(
                "data/old/kext_albedo_WD_MW_5.5A_30_D03.all_modified",
                format="ascii.commented_header",
            )
            ax.plot(
                b["lambda"],
                b["C_ext/H"] / 6.622e-22,
                "k:",
                label="MW R(V)=5.5; sizedist=A (Draine 2003)",
                alpha=0.5,
            )

            b = Table.read(
                "data/old/kext_albedo_WD_MW_5.5B_30.dat_modified",
                format="ascii.commented_header",
            )
            ax.plot(
                b["lambda"],
                b["C_ext/H"] / 4.789e-22,
                "k-.",
                label="MW R(V)=5.5; sizedist=B (Weingartner & Draine 2001)",
                alpha=0.5,
            )

        # get the average extinction curve
        if args.ave:
            ave_extdata = AverageExtData(extdatas, alav=True)
            ave_extdata.plot(
                ax,
                color="k",
                fontsize=fontsize,
                legend_key="IRS",
                rebin_fac=args.rebin_fac,
                legend_label="Average (this work)",
            )
            ave_extdata.save(args.filelist.replace(".dat", "_ave.fits"))

        # ax.legend(fontsize=fontsize)
    else:
        ax.set_xlim(1.0, 100.0)
        ax.set_ylim(-6, -1.0)
        ax.set_ylabel(r"$E(\lambda - V)$", fontsize=1.3 * fontsize)
        # ax.legend(fontsize=12, ncol=4)

    ax.set_xlabel(r"$\lambda$ [$\mu m$]")

    ax.tick_params("both", length=10, width=2, which="major")
    ax.tick_params("both", length=5, width=1, which="minor")

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    ax.xaxis.set_minor_formatter(ScalarFormatter())
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.set_xticks([1.0, 10.0])
    ax.set_xticks([2, 3, 4, 5, 6, 7, 8, 15.0, 20.0, 30.0, 40.0], minor=True)

    fig.tight_layout()

    save_str = "_mext"
    if args.png:
        fig.savefig(args.filelist.replace(".dat", save_str + ".png"))
    elif args.pdf:
        fig.savefig(args.filelist.replace(".dat", save_str + ".pdf"))
    else:
        pyplot.show()
