# plot to plot the silicate strength versus other properties of the sightlines

import argparse
import math
import os.path

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib
import astropy.units as u
from astropy.table import Table

from measure_extinction.extdata import ExtData

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filelist", help="file with list of curves to plot")
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    filename = args.filelist

    # data
    f = open(filename, "r")
    file_lines = list(f)
    extnames = []
    extdatas = []
    extdatas_fm90 = []
    for line in file_lines:
        if (line.find("#") != 0) & (len(line) > 0):
            name = line.rstrip()
            extnames.append(name.split("_")[0])
            bfile = f"fits/{name}"
            text = ExtData(filename=bfile)
            extdatas.append(text)

            uvfname = bfile.replace(".fits", "_FM90.fits")
            if os.path.isfile(uvfname):
                text = ExtData(filename=uvfname)
            else:
                text = None
            extdatas_fm90.append(text)

    # plots
    fontsize = 14

    font = {"size": fontsize}

    matplotlib.rc("font", **font)

    matplotlib.rc("lines", linewidth=1)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("xtick.minor", width=2)
    matplotlib.rc("ytick.major", width=2)
    matplotlib.rc("ytick.minor", width=2)

    figsize = (9, 9)
    fig, fax = pyplot.subplots(nrows=2, ncols=2, figsize=figsize)
    ax = [fax[0, 0], fax[0, 1], fax[1, 0], fax[1, 1]]

    n_ext = len(extdatas)
    sil_amp = np.full((n_ext), 0.0)
    sil_amp_unc = np.full((n_ext), 0.0)
    sil_width = np.full((n_ext), 0.0)
    sil_width_unc = np.full((n_ext), 0.0)
    sil_area = np.full((n_ext), 0.0)
    sil_area_unc = np.full((n_ext), 0.0)
    nuv_amp = np.full((n_ext), 0.0)
    nuv_amp_unc = np.full((n_ext), 0.0)
    nuv_area = np.full((n_ext), 0.0)
    nuv_area_unc = np.full((n_ext), 0.0)
    avs = np.full((n_ext), 0.0)
    avs_unc = np.full((n_ext), 0.0)
    ebvs = np.full((n_ext), 0.0)
    ebvs_unc = np.full((n_ext), 0.0)
    for k, cext in enumerate(extdatas):
        avs[k] = cext.columns["AV"][0]
        avs_unc[k] = cext.columns["AV"][1]

        (indxs,) = np.where(
            (cext.waves["BAND"] > 0.4 * u.micron)
            & (cext.waves["BAND"] < 0.5 * u.micron)
        )
        ebvs[k] = cext.exts["BAND"][indxs[0]]
        ebvs_unc[k] = cext.uncs["BAND"][indxs[0]]

        samp = cext.p92_p50_fit["SIL1_AMP"]
        swidth = cext.p92_p50_fit["SIL1_WIDTH"]
        sil_amp[k] = samp[0]
        sil_amp_unc[k] = samp[1]
        sil_width[k] = swidth[0]
        sil_width_unc[k] = swidth[1]
        sil_area[k] = math.pi * samp[0] / (2.0 * swidth[0])
        sil_area_unc[k] = (0.5 * (samp[1] + samp[2]) / samp[0]) ** 2 + (
            0.5 * (swidth[1] + swidth[2]) / samp[0]
        ) ** 2

        cext_fm90 = extdatas_fm90[k]
        if cext_fm90 is not None:
            samp = cext_fm90.fm90_p50_fit["C3"]
            swidth = cext_fm90.fm90_p50_fit["GAMMA"]
            nuv_amp[k] = samp[0]
            nuv_amp_unc[k] = samp[1]
            nuv_area[k] = math.pi * samp[0] / (2.0 * swidth[0])
            nuv_area_unc[k] = (0.5 * (samp[1] + samp[2]) / samp[0]) ** 2 + (
                0.5 * (swidth[1] + swidth[2]) / samp[0]
            ) ** 2

    sil_area_unc = np.sqrt(sil_area_unc) * sil_area
    nuv_area_unc = np.sqrt(nuv_area_unc) * nuv_area

    rvs = avs / ebvs
    rvs_unc = (avs_unc / avs) ** 2 + (ebvs_unc / ebvs) ** 2
    rvs_unc = rvs * np.sqrt(rvs_unc)

    # output some info
    a = Table()
    a["name"] = extnames
    a["AV"] = avs
    a["RV"] = rvs
    a["Sil_amp"] = sil_amp
    a["Sil_width"] = sil_width
    a["Sil_area"] = sil_area
    print(a)

    # uncs are way over estimated as the amplitude and width are very well correlated
    # update once the samples are available
    # print(sil_amp_unc)
    sil_area_unc = np.full((n_ext), 0.0)

    gindxs = sil_area < 0.01
    bindxs = sil_area > 0.01

    # R(V) versus A(V)
    ax[0].errorbar(
        rvs[gindxs], avs[gindxs], xerr=rvs_unc[gindxs], yerr=avs_unc[gindxs], fmt="ko"
    )
    ax[0].errorbar(
        rvs[bindxs],
        avs[bindxs],
        xerr=rvs_unc[bindxs],
        yerr=avs_unc[bindxs],
        fmt="ko",
        markerfacecolor="none",
    )
    ax[0].set_xlabel(r"R(V) = A(V)/E(B-V)")
    ax[0].set_ylabel(r"A(V)")
    ax[0].tick_params("both", length=10, width=2, which="major")
    ax[0].tick_params("both", length=5, width=1, which="minor")

    # R(V) versus silicate
    ax[1].errorbar(
        rvs[gindxs],
        sil_area[gindxs],
        xerr=rvs_unc[gindxs],
        yerr=sil_area_unc[gindxs],
        fmt="ko",
    )
    ax[1].set_xlabel(r"R(V) = A(V)/E(B-V)")
    ax[1].set_ylabel(r"Silicate 10 $\mu$m area")
    ax[1].tick_params("both", length=10, width=2, which="major")
    ax[1].tick_params("both", length=5, width=1, which="minor")

    # A(V) versus A(sil)/A(V)
    ax[3].errorbar(
        avs[gindxs],
        sil_amp[gindxs],
        xerr=avs_unc[gindxs],
        yerr=sil_amp_unc[gindxs],
        fmt="ko",
    )
    ax[3].set_xlabel(r"A(V)")
    ax[3].set_ylabel(r"$A(sil)/A(V) \approx \tau(sil)/A(V)$")
    ax[3].tick_params("both", length=10, width=2, which="major")
    ax[3].tick_params("both", length=5, width=1, which="minor")

    # silicate verus 2175
    uvindxs = (sil_area < 0.01) & (nuv_area > 0.0)
    ax[2].errorbar(
        sil_area[uvindxs],
        nuv_area[uvindxs],
        xerr=sil_area_unc[uvindxs],
        yerr=nuv_area_unc[uvindxs],
        fmt="ko",
    )
    ax[2].set_yscale("linear")
    ax[2].set_xscale("linear")
    ax[2].set_xlabel(r"Silicate 10 $\mu$m area")
    ax[2].set_ylabel(r"2175 $\AA$ bump area")
    # ax[2].set_xlim(kxrange)

    ax[2].tick_params("both", length=10, width=2, which="major")
    ax[2].tick_params("both", length=5, width=1, which="minor")

    fig.tight_layout()

    save_str = "_silicate"
    if args.png:
        fig.savefig(args.filelist.replace(".dat", save_str + ".png"))
    elif args.pdf:
        fig.savefig(args.filelist.replace(".dat", save_str + ".pdf"))
    else:
        pyplot.show()
