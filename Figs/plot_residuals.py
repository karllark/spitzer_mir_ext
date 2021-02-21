import argparse

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib
from matplotlib.ticker import ScalarFormatter

import astropy.units as u

from measure_extinction.merge_obsspec import _wavegrid
from measure_extinction.extdata import ExtData

from utils.G21 import G21_drude_asym as G21

if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    fontsize = 14

    font = {"size": fontsize}

    matplotlib.rc("font", **font)

    matplotlib.rc("lines", linewidth=1)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("xtick.minor", width=2)
    matplotlib.rc("ytick.major", width=2)
    matplotlib.rc("ytick.minor", width=2)

    figsize = (12, 4)
    fig, ax = pyplot.subplots(figsize=figsize)

    # New measurements
    avefilenames = [
        "data/all_ext_14oct20_diffuse_ave_POWLAW2DRUDE.fits",
        # "fits_good_18aug20/hd283809_hd064802_ext_POWLAW2DRUDE.fits",
        # "fits_good_18aug20/hd029647_hd195986_ext_POWLAW2DRUDE.fits",
    ]
    pcol = ["b", "g", "c"]
    psym = ["o", "s", "^"]
    pline = ["-", ":", "-."]
    plabel = ["diffuse", "HD283809", "HD029647"]

    for i, avefilename in enumerate(avefilenames):
        # IR Fit
        obsext = ExtData()
        obsext.read(avefilename)

        if obsext.type == "elx":
            obsext.trans_elv_alav()

        obsext_wave = obsext.waves["BAND"].value

        obsext_ext = obsext.exts["BAND"]
        obsext_ext_uncs = obsext.uncs["BAND"]

        gindxs_IRS = np.where(obsext.npts["IRS"] > 0)
        obsext_IRS_wave = obsext.waves["IRS"][gindxs_IRS].value
        obsext_IRS_ext = obsext.exts["IRS"][gindxs_IRS]
        obsext_IRS_uncs = obsext.uncs["IRS"][gindxs_IRS]

        G21_p50 = G21(
            scale=obsext.g21_p50_fit["SCALE"][0],
            alpha=obsext.g21_p50_fit["ALPHA"][0],
            sil1_amp=obsext.g21_p50_fit["SIL1_AMP"][0],
            sil1_center=obsext.g21_p50_fit["SIL1_CENTER"][0],
            sil1_fwhm=obsext.g21_p50_fit["SIL1_FWHM"][0],
            sil1_asym=obsext.g21_p50_fit["SIL1_ASYM"][0],
            sil2_amp=obsext.g21_p50_fit["SIL2_AMP"][0],
            sil2_center=obsext.g21_p50_fit["SIL2_CENTER"][0],
            sil2_fwhm=obsext.g21_p50_fit["SIL2_FWHM"][0],
            sil2_asym=obsext.g21_p50_fit["SIL2_ASYM"][0],
        )

        # gvals = obsext_wave > 1.0
        # ax.errorbar(
        #    obsext_wave[gvals],
        #    obsext_ext[gvals] - G21_p50(obsext_wave[gvals] * u.micron),
        #    yerr=obsext_ext_uncs[gvals],
        #    fmt=pcol[i] + psym[i],
        #    markersize=10,
        #    markeredgewidth=1.0,
        #    alpha=0.5,
        # )

        # IRS
        IRS_resid = obsext_IRS_ext - G21_p50(obsext_IRS_wave * u.micron)
        IRS_wave = obsext_IRS_wave
        ax.plot(
            IRS_wave,
            IRS_resid,
            pcol[i] + pline[i],
            lw=2,
            alpha=0.65,
            label="Diffuse",
        )

        # rebin IRS
        wrange = [5.0, 36.0]
        res = 25
        full_wave, full_wave_min, full_wave_max = _wavegrid(res, wrange)
        n_waves = len(full_wave)
        full_flux = np.zeros((n_waves), dtype=float)
        full_unc = np.zeros((n_waves), dtype=float)
        full_npts = np.zeros((n_waves), dtype=int)

        cwaves = obsext_IRS_wave
        cfluxes = obsext_IRS_ext
        cuncs = obsext_IRS_uncs
        for k in range(n_waves):
            (indxs,) = np.where(
                (cwaves >= full_wave_min[k]) & (cwaves < full_wave_max[k])
            )
            if len(indxs) > 0:
                # weights = 1.0 / np.square(cuncs[indxs])
                weights = 1.0
                full_flux[k] += np.sum(weights * cfluxes[indxs])
                full_unc[k] += np.sum(1.0 / np.square(cuncs[indxs]))
                full_npts[k] += len(indxs)

        findxs = full_npts > 0
        full_flux[findxs] /= full_npts[findxs]
        full_unc[findxs] = np.sqrt(1.0 / full_unc[findxs])

        ax.errorbar(
            full_wave[findxs],
            full_flux[findxs] - G21_p50(full_wave[findxs] * u.micron),
            yerr=full_unc[findxs],
            fmt="bo",  # pcol[i] + psym[i],
            markersize=5,
            markeredgewidth=1.0,
            alpha=1.0,
            label="Diffuse R=25",
        )

        # average residuals
        minwaves = [5.0, 10.0, 20.0]
        maxwaves = [10.0, 20.0, 40.0]
        for minw, maxw in zip(minwaves, maxwaves):
            totstd = np.std(IRS_resid[(IRS_wave >= minw) & (IRS_wave < maxw)])
            avewave = 0.5 * (minw + maxw)
            ax.text(
                avewave,
                -0.025,
                rf"$\sigma = {totstd:.4f}$",
                ha="center",
                bbox=dict(facecolor="white", edgecolor="white", alpha=0.75),
            )
            ax.plot([minw, maxw, maxw], [-0.02, -0.02, -0.023], "k-", alpha=0.5, lw=2)

    ax.set_yscale("linear")
    ax.set_xlabel(r"$\lambda$ [$\mu m$]")
    ax.tick_params("both", length=10, width=2, which="major")
    ax.tick_params("both", length=5, width=1, which="minor")
    ax.set_ylabel(r"$A(\lambda)/A(V)$", fontsize=1.3 * fontsize)

    ax.plot([5.0, 40.0], [0.0, 0.0], "k--", alpha=0.75, lw=2)
    ax.plot([5.0, 40.0], [0.002, 0.002], "k:", alpha=0.5, lw=2)
    linewaves = [6.2, 6.85, 7.7, 11.1]
    for cwave in linewaves:
        ax.plot([cwave, cwave], [0.0, 0.01], "k-.", alpha=0.5, lw=2)
        ax.text(
            cwave,
            0.011,
            rf"{cwave} $\mu$m",
            ha="center",
            rotation="vertical",
            fontsize=0.8 * fontsize,
            alpha=0.7,
        )

    # ax.annotate("6.2", )
    ax.set_xscale("log")
    ax.set_xlim(5.0, 40.0)
    ax.set_ylim(-0.03, 0.03)
    # ax.yaxis.tick_right()
    # ax.yaxis.set_label_position("right")
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(ScalarFormatter())
    ax.yaxis.set_major_formatter(ScalarFormatter())

    ax.legend(loc="upper right")

    fig.tight_layout()

    save_fname = "irext_residuals"
    if args.png:
        fig.savefig(save_fname + ".png")
    elif args.pdf:
        fig.savefig(save_fname + ".pdf")
    else:
        pyplot.show()
