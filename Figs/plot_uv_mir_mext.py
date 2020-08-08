import argparse
import os.path

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib

# from astropy.table import Table
import astropy.units as u

# from calc_ext import P92_Elv
from dust_extinction.shapes import FM90
from utils.P92_mod import P92_mod as P92
from measure_extinction.extdata import ExtData
from utils.G21 import G20_drude_asym as G21


def plot_all_ext(
    ax, kxrange, kyrange, normvals=None, yoffset_factor=0.0, annotate_key=None
):
    """
    plot all the extintion info on the specified plot
    """
    # sindxs = np.argsort(avs)
    sindxs = np.arange(len(avs))

    ann_wave_range = [5.0, 10.0] * u.micron
    col_vals = ["b", "g"]  # , "r", "m", "c", "y"]
    lin_vals = ["--", ":", "-."]
    n_cols = len(col_vals)

    mod_x = np.logspace(0.0, 2.0, 200) * u.micron
    mod_x_g21 = np.logspace(0.1, np.log10(39.), 200) * u.micron
    mod_x_fm90 = np.logspace(-1.0, -0.5, 200) * u.micron
    for i in range(len(extnames)):
        k = sindxs[i]

        if normvals is not None:
            normval = normvals[k]
        else:
            normval = 1.0

        # plot the extinction curves
        if not args.modonly:
            extdatas[k].plot(
                ax,
                color=col_vals[i % n_cols],
                alax=True,
                normval=normval,
                yoffset=i * yoffset_factor,
                alpha=1.0,
                rebin_fac=args.rebin_fac,
                fontsize=fontsize,
                annotate_key=annotate_key,
                annotate_yoffset=0.025,
                annotate_text=extnames[k].split("_")[0],
                annotate_wave_range=ann_wave_range,
            )

        if args.models:

            # plot the best fit P92 model
            if hasattr(extdatas[k], 'p92_p50_fit'):
                P92_best = P92(
                    BKG_amp=extdatas[k].p92_p50_fit["BKG_AMP"][0],
                    BKG_lambda=extdatas[k].p92_p50_fit["BKG_LAMBDA"][0],
                    BKG_width=extdatas[k].p92_p50_fit["BKG_WIDTH"][0],
                    FUV_amp=extdatas[k].p92_p50_fit["FUV_AMP"][0],
                    FUV_lambda=extdatas[k].p92_p50_fit["FUV_LAMBDA"][0],
                    FUV_n=extdatas[k].p92_p50_fit["FUV_N"][0],
                    FUV_b=extdatas[k].p92_p50_fit["FUV_B"][0],
                    NUV_amp=extdatas[k].p92_p50_fit["NUV_AMP"][0],
                    NUV_lambda=extdatas[k].p92_p50_fit["NUV_LAMBDA"][0],
                    NUV_width=extdatas[k].p92_p50_fit["NUV_WIDTH"][0],
                    SIL1_amp=extdatas[k].p92_p50_fit["SIL1_AMP"][0],
                    SIL1_lambda=extdatas[k].p92_p50_fit["SIL1_LAMBDA"][0],
                    SIL1_width=extdatas[k].p92_p50_fit["SIL1_WIDTH"][0],
                    SIL2_amp=extdatas[k].p92_p50_fit["SIL2_AMP"][0],
                    SIL2_lambda=extdatas[k].p92_p50_fit["SIL2_LAMBDA"][0],
                    SIL2_width=extdatas[k].p92_p50_fit["SIL2_WIDTH"][0],
                    FIR_amp=extdatas[k].p92_p50_fit["FIR_AMP"][0],
                    FIR_lambda=extdatas[k].p92_p50_fit["FIR_LAMBDA"][0],
                    FIR_width=extdatas[k].p92_p50_fit["FIR_WIDTH"][0],
                )

                ax.plot(
                    mod_x,
                    P92_best(mod_x) / normval + i * yoffset_factor,
                    lin_vals[i % 3],
                    color=col_vals[i % n_cols],
                    alpha=0.5,
                )

            if hasattr(extdatas[k], 'g21_best_fit'):
                # best fit G21 model
                if extdatas[k] is not None:
                    G21_best = G21(
                        scale=extdatas[k].g21_best_fit["SCALE"],
                        alpha=extdatas[k].g21_best_fit["ALPHA"],
                        sil1_amp=extdatas[k].g21_best_fit["SIL1_AMP"],
                        sil1_center=extdatas[k].g21_best_fit["SIL1_CENTER"],
                        sil1_fwhm=extdatas[k].g21_best_fit["SIL1_FWHM"],
                        sil1_asym=extdatas[k].g21_best_fit["SIL1_ASYM"],
                        # sil2_amp=0.0,
                        sil2_amp=extdatas[k].g21_best_fit["SIL2_AMP"],
                        sil2_center=extdatas[k].g21_best_fit["SIL2_CENTER"],
                        sil2_fwhm=extdatas[k].g21_best_fit["SIL2_FWHM"],
                        sil2_asym=extdatas[k].g21_best_fit["SIL2_ASYM"],
                    )

                ax.plot(
                    mod_x_g21,
                    G21_best(mod_x_g21) / normval + i * yoffset_factor,
                    lin_vals[i % 3],
                    color=col_vals[i % n_cols],
                    alpha=0.5,
                )

            if extdatas_fm90[k] is not None:
                if hasattr(extdatas_fm90[k], 'fm90_best_fit'):
                    # best fit FM90 model
                    if extdatas_fm90[k] is not None:
                        FM90_best = FM90(
                            C1=extdatas_fm90[k].fm90_best_fit["C1"],
                            C2=extdatas_fm90[k].fm90_best_fit["C2"],
                            C3=extdatas_fm90[k].fm90_best_fit["C3"],
                            C4=extdatas_fm90[k].fm90_best_fit["C4"],
                            xo=extdatas_fm90[k].fm90_best_fit["XO"],
                            gamma=extdatas_fm90[k].fm90_best_fit["GAMMA"],
                        )

                        ax.plot(
                            mod_x_fm90,
                            FM90_best(mod_x_fm90) / normval + i * yoffset_factor,
                            lin_vals[i % 3],
                            color=col_vals[i % n_cols],
                            alpha=0.5,
                        )

    ax.set_yscale("linear")
    ax.set_xscale("linear")
    ax.set_xlim(kxrange)
    ax.set_ylabel(r"$A(\lambda)/A(V)$", fontsize=1.3 * fontsize)

    ax.set_xlabel(r"$\lambda$ [$\mu m$]")

    ax.tick_params("both", length=10, width=2, which="major")
    ax.tick_params("both", length=5, width=1, which="minor")


if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("filelist", help="file with list of curves to plot")
    parser.add_argument(
        "--rebin_fac", type=int, default=None, help="rebin factor for spectra"
    )
    parser.add_argument(
        "--models", help="plot the best fit models", action="store_true"
    )
    parser.add_argument(
        "--modonly", help="only plot the best fit models", action="store_true"
    )
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    filename = args.filelist

    f = open(filename, "r")
    file_lines = list(f)
    extnames = []
    extdatas = []
    extdatas_fm90 = []
    avs = []

    normtype = "IUE"
    norm_wave_range = [0.25, 0.30] * u.micron
    normvals = []

    for line in file_lines:
        if (line.find("#") != 0) & (len(line) > 0):
            name = line.rstrip()
            extnames.append(name)
            bfilename = f"fits/{name}"
            text = ExtData(filename=bfilename)
            extdatas.append(text)
            avs.append(text.columns["AV"][0])

            # determine the extinction in the near-UV
            # useful for sorting to make a pretty plot
            if "IUE" in text.exts.keys():
                (gindxs,) = np.where(
                    (text.npts[normtype] > 0)
                    & (
                        (text.waves[normtype] >= norm_wave_range[0])
                        & (text.waves[normtype] <= norm_wave_range[1])
                    )
                )
                normvals.append(
                    np.average(
                        (text.exts[normtype][gindxs]) / float(text.columns["AV"][0])
                        + 1.0
                    )
                )
            else:
                normvals.append(1.0)
        else:
            normvals.append(1.0)

    # print(normvals)
    # normvals = np.array(normvals)
    # sindxs = np.flip(np.argsort(normvals))
    # normvals = normvals[sindxs]
    # extnames = np.array(extnames)[sindxs]

    extdatas = []
    avs = []

    for extname in extnames:
        bfilename = f"fits/{extname}"
        text = ExtData(filename=bfilename)
        extdatas.append(text)
        avs.append(text.columns["AV"][0])

        fm90_filename = bfilename.replace(".fits", "_FM90.fits")
        if os.path.isfile(fm90_filename):
            textfm90 = ExtData(filename=fm90_filename)
            extdatas_fm90.append(textfm90)
        else:
            extdatas_fm90.append(None)

    fontsize = 18

    font = {"size": fontsize}

    matplotlib.rc("font", **font)

    matplotlib.rc("lines", linewidth=1)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("xtick.minor", width=2)
    matplotlib.rc("ytick.major", width=2)
    matplotlib.rc("ytick.minor", width=2)

    figsize = (12, 10)
    fig, ax = pyplot.subplots(nrows=1, ncols=2, figsize=figsize)

    plot_all_ext(
        ax[1],
        kxrange=[1.0, 40.0],
        kyrange=[-6.0, -0.5],
        yoffset_factor=0.1,
        annotate_key="IRS",
    )
    plot_all_ext(
        ax[0],
        kxrange=[0.115, 0.33],
        kyrange=[1.0, 10.0],
        normvals=normvals,
        yoffset_factor=0.5,
    )

    ax[1].set_ylim(-0.1, 1.75)
    ax[0].set_ylim(0.0, 9.5)
    ax[1].set_ylabel(r"$A(\lambda)/A(V) + constant$")
    ax[0].set_ylabel(r"$A(\lambda)/A(0.27~\mu m) + constant$")

    ax[1].yaxis.set_label_position("right")
    ax[1].yaxis.tick_right()

    fig.tight_layout()

    save_str = "_mext_uv_mir"
    if args.png:
        fig.savefig(args.filelist.replace(".dat", save_str + ".png"))
    elif args.pdf:
        fig.savefig(args.filelist.replace(".dat", save_str + ".pdf"))
    else:
        pyplot.show()
