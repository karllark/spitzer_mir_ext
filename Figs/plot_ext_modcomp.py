import argparse

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib
from matplotlib.ticker import ScalarFormatter
from matplotlib.lines import Line2D

import astropy.units as u

from dust_extinction.grain_models import D03, ZDA04, J13
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

    figsize = (8, 6)
    fig, ax = pyplot.subplots(figsize=figsize)

    # New measurements
    avefilename = "data/all_ext_14oct20_diffuse_ave_POWLAW2DRUDE.fits"

    # IR Fit
    obsext = ExtData()
    obsext.read(avefilename)

    obsext_wave = obsext.waves["BAND"].value

    obsext_ext = obsext.exts["BAND"]
    obsext_ext_uncs = obsext.uncs["BAND"]

    gindxs_IRS = np.where(obsext.npts["IRS"] > 0)
    obsext_IRS_wave = obsext.waves["IRS"][gindxs_IRS].value
    obsext_IRS_ext = obsext.exts["IRS"][gindxs_IRS]
    obsext_IRS_uncs = obsext.uncs["IRS"][gindxs_IRS]

    ax.errorbar(
        obsext_wave,
        obsext_ext,
        yerr=obsext_ext_uncs,
        fmt="bo",
        markersize=10,
        markeredgewidth=1.0,
        alpha=0.5,
    )

    # IRS
    # ax.plot(
    #     obsext_IRS_wave, obsext_IRS_ext, "b-", lw=2, alpha=0.65,
    # )

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
        (indxs,) = np.where((cwaves >= full_wave_min[k]) & (cwaves < full_wave_max[k]))
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
        full_flux[findxs],
        yerr=full_unc[findxs],
        fmt="bo",  # pcol[i] + psym[i],
        markersize=5,
        markeredgewidth=1.0,
        alpha=1.0,
    )

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

    mod_x = np.logspace(np.log10(1.0), np.log10(39.0), num=1000) * u.micron
    ax.plot(mod_x, G21_p50(mod_x), "k-", lw=2, alpha=0.65, label="G21 Fit")

    dmod = D03()
    ax.plot(mod_x, dmod(mod_x), "g:", lw=2, label="D03 MWRV31")

    zmod = ZDA04()
    ax.plot(mod_x, zmod(mod_x), "m--", lw=2, label="ZDA04 MWRV31")

    jmod = J13()
    ax.plot(mod_x, jmod(mod_x), "c-.", lw=2, label="J11 MWRV31")

    ax.set_yscale("linear")
    ax.set_xlabel(r"$\lambda$ [$\mu m$]")
    ax.tick_params("both", length=10, width=2, which="major")
    ax.tick_params("both", length=5, width=1, which="minor")
    ax.set_ylabel(r"$A(\lambda)/A(V)$", fontsize=1.3 * fontsize)

    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlim(1.0, 40.0)
    ax.set_ylim(0.007, 0.4)
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_major_formatter(ScalarFormatter())

    # custom legend
    custom_lines = [
        Line2D(
            [0],
            [0],
            color="w",
            marker="o",
            markersize=5,
            markerfacecolor="b",
            alpha=0.5,
        ),
        Line2D(
            [0],
            [0],
            color="w",
            marker="o",
            markerfacecolor="b",
            markersize=5,
            alpha=1.0,
        ),
        Line2D([0], [0], color="k", lw=2, alpha=0.65),
        Line2D([0], [0], color="g", lw=2, linestyle=":"),
        Line2D([0], [0], color="m", lw=2, linestyle="--"),
        Line2D([0], [0], color="c", lw=2, linestyle="-."),
    ]
    ax.legend(  # w/o "Diffuse R=25",
        custom_lines, ["Diffuse phot", "Diffuse R=25", "G21 Fit", "D03", "ZDA04", "J13"],
    )

    ax.xaxis.set_minor_formatter(ScalarFormatter())
    ax.set_xticks([2, 3, 4, 5, 6, 7, 8, 15.0, 20.0, 30.0], minor=True)

    fig.tight_layout()

    save_fname = "diffuse_ext_modcomp"
    if args.png:
        fig.savefig(save_fname + ".png")
    elif args.pdf:
        fig.savefig(save_fname + ".pdf")
    else:
        pyplot.show()
