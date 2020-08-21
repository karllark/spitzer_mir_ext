import argparse

import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib
from matplotlib.ticker import ScalarFormatter
from matplotlib.lines import Line2D

import astropy.units as u
from astropy.table import QTable

from dust_extinction.parameter_averages import F19
from dust_extinction.shapes import FM90
from measure_extinction.merge_obsspec import _wavegrid
from measure_extinction.extdata import ExtData

from utils.G21 import G20_drude_asym as G21

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

    figsize = (14, 6)
    fig, ax = pyplot.subplots(nrows=1, ncols=2, figsize=figsize)

    # F19 R(V) average model
    mod_x = np.logspace(np.log10(0.115), np.log10(2.5), num=1000) * u.micron
    F19_Rv = F19(Rv=3.1)
    for cax in ax:
        cax.plot(mod_x, F19_Rv(mod_x), "k:", lw=2, alpha=0.65, label="F19 R(V) = 3.1")

    # New measurements
    avefilenames = [
        "data/all_ext_18feb20_diffuse_ave_POWLAW2DRUDE.fits",
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

        # UV fit
        obsext2 = ExtData()
        if plabel[i] == "diffuse":
            obsext2.read(avefilename.replace("POWLAW2DRUDE", "FM90"))
        else:
            obsext2.read(avefilename.replace(".fits", "_FM90.fits"))

        obsext_wave = obsext.waves["BAND"].value

        obsext_ext = obsext.exts["BAND"]
        obsext_ext_uncs = obsext.uncs["BAND"]

        gindxs_IRS = np.where(obsext.npts["IRS"] > 0)
        obsext_IRS_wave = obsext.waves["IRS"][gindxs_IRS].value
        obsext_IRS_ext = obsext.exts["IRS"][gindxs_IRS]
        obsext_IRS_uncs = obsext.uncs["IRS"][gindxs_IRS]

        gindxs_IUE = np.where(obsext2.npts["IUE"] > 0)
        obsext_IUE_wave = obsext2.waves["IUE"][gindxs_IUE].value
        obsext_IUE_ext = obsext2.exts["IUE"][gindxs_IUE]
        obsext_IUE_uncs = obsext2.uncs["IUE"][gindxs_IUE]

        ax[1].errorbar(
            obsext_wave,
            obsext_ext,
            yerr=obsext_ext_uncs,
            fmt=pcol[i] + psym[i],
            markersize=10,
            markeredgewidth=1.0,
            alpha=0.5,
        )

        # IRS
        ax[1].plot(
            obsext_IRS_wave, obsext_IRS_ext, pcol[i] + pline[i], lw=2, alpha=0.65,
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

        ax[1].errorbar(
            full_wave[findxs],
            full_flux[findxs],
            yerr=full_unc[findxs],
            fmt="bo",  # pcol[i] + psym[i],
            markersize=5,
            markeredgewidth=1.0,
            alpha=1.0,
        )

        # Opt photometry
        ax[0].errorbar(
            obsext_wave,
            obsext_ext,
            yerr=obsext_ext_uncs,
            fmt=pcol[i] + psym[i],
            markersize=10,
            markeredgewidth=1.0,
            alpha=0.5,
        )

        # IUE
        ax[0].plot(
            obsext_IUE_wave, obsext_IUE_ext, pcol[i] + pline[i], lw=2, alpha=0.85,
        )

        G21_best = G21(
            scale=obsext.g21_best_fit["SCALE"],
            alpha=obsext.g21_best_fit["ALPHA"],
            sil1_amp=obsext.g21_best_fit["SIL1_AMP"],
            sil1_center=obsext.g21_best_fit["SIL1_CENTER"],
            sil1_fwhm=obsext.g21_best_fit["SIL1_FWHM"],
            sil1_asym=obsext.g21_best_fit["SIL1_ASYM"],
            sil2_amp=obsext.g21_best_fit["SIL2_AMP"],
            sil2_center=obsext.g21_best_fit["SIL2_CENTER"],
            sil2_fwhm=obsext.g21_best_fit["SIL2_FWHM"],
            sil2_asym=obsext.g21_best_fit["SIL2_ASYM"],
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
        ax[1].plot(
            mod_x, G21_p50(mod_x), "k" + pline[i], lw=2, alpha=0.65, label="G21 Fit"
        )

        g21_comps = G21_p50.copy()
        g21_comps.sil1_amp = 0.0
        ax[1].plot(mod_x, g21_comps(mod_x), "k--", alpha=0.5)

        g21_comps = G21_p50.copy()
        g21_comps.sil2_amp = 0.0
        ax[1].plot(mod_x, g21_comps(mod_x), "k--", alpha=0.5)

        g21_comps = G21_p50.copy()
        g21_comps.sil1_amp = 0.0
        g21_comps.sil2_amp = 0.0
        ax[1].plot(mod_x, g21_comps(mod_x), "k--", alpha=0.5)

        # write ext to table
        gphot = obsext_wave > 1.0
        rebintab = QTable()
        twaves = np.concatenate([obsext_wave[gphot], full_wave[findxs]]) * u.micron
        text = np.concatenate([obsext_ext[gphot], full_flux[findxs]])
        tunc = np.concatenate([obsext_ext_uncs[gphot], full_unc[findxs]])
        nt = len(twaves)
        rebintab["wave1"] = twaves[0 : nt // 2]
        rebintab["ext1"] = text[0 : nt // 2]
        rebintab["unc1"] = tunc[0 : nt // 2]
        rebintab["fit1"] = G21_p50(rebintab["wave1"])
        rebintab["wave2"] = twaves[nt // 2 : nt]
        rebintab["ext2"] = text[nt // 2 : nt]
        rebintab["unc2"] = tunc[nt // 2 : nt]
        rebintab["fit2"] = G21_p50(rebintab["wave2"])
        rebintab.write(
            "test.tex",
            formats={
                "wave1": "%4.2f",
                "ext1": "%4.4f",
                "unc1": "%4.4f",
                "fit1": "%4.4f",
                "wave2": "%4.2f",
                "ext2": "%4.4f",
                "unc2": "%4.4f",
                "fit2": "%4.4f",
            },
            overwrite=True,
        )

        # UV
        FM90_p50 = FM90(
            C1=obsext2.fm90_p50_fit["C1"][0],
            C2=obsext2.fm90_p50_fit["C2"][0],
            C3=obsext2.fm90_p50_fit["C3"][0],
            C4=obsext2.fm90_p50_fit["C4"][0],
            xo=obsext2.fm90_p50_fit["XO"][0],
            gamma=obsext2.fm90_p50_fit["GAMMA"][0],
        )

        mod_x = np.logspace(np.log10(0.1), np.log10(0.3), num=1000) * u.micron

        ax[0].plot(
            mod_x, FM90_p50(mod_x), "k" + pline[i], lw=2, alpha=0.65, label="FM90 Fit"
        )

        fm90_comps = FM90_p50.copy()
        fm90_comps.C3 = 0.0
        ax[0].plot(mod_x, fm90_comps(mod_x), "k--", alpha=0.5)

        fm90_comps = FM90_p50.copy()
        fm90_comps.C4 = 0.0
        ax[0].plot(mod_x, fm90_comps(mod_x), "k--", alpha=0.5)

        fm90_comps = FM90_p50.copy()
        fm90_comps.C3 = 0.0
        fm90_comps.C4 = 0.0
        ax[0].plot(mod_x, fm90_comps(mod_x), "k--", alpha=0.5)

    for i in range(2):
        ax[i].set_yscale("linear")
        ax[i].set_xlabel(r"$\lambda$ [$\mu m$]")
        ax[i].tick_params("both", length=10, width=2, which="major")
        ax[i].tick_params("both", length=5, width=1, which="minor")
        ax[i].set_ylabel(r"$A(\lambda)/A(V)$", fontsize=1.3 * fontsize)

        # ax[i].legend(fontsize=0.75 * fontsize)

    # finishing plot details
    ax[0].set_xlim(0.1, 0.6)
    ax[0].set_xscale("log")
    ax[0].set_ylim(0.0, 4.0)
    ax[0].xaxis.set_major_formatter(ScalarFormatter())
    ax[0].xaxis.set_minor_formatter(ScalarFormatter())

    ax[1].set_yscale("log")
    ax[1].set_xscale("log")
    ax[1].set_xlim(1.0, 40.0)
    ax[1].set_ylim(0.01, 0.4)
    ax[1].yaxis.tick_right()
    ax[1].yaxis.set_label_position("right")
    ax[1].xaxis.set_major_formatter(ScalarFormatter())
    ax[1].yaxis.set_major_formatter(ScalarFormatter())

    # custom legend
    # custom legend
    custom_lines = [
        Line2D([0], [0], color="b", marker="o", markersize=10, alpha=0.5),
        Line2D([0], [0], color="k", lw=2, alpha=0.65),
        Line2D([0], [0], color="k", lw=2, alpha=0.5, linestyle="--"),
        Line2D([0], [0], color="k", lw=2, linestyle=":", alpha=0.65),
    ]
    ax[0].legend(
        custom_lines, ["Diffuse", "FM90 Fit", "FM90 Components", "F19 R(V)=3.1"]
    )

    # custom legend
    custom_lines = [
        Line2D([0], [0], color="b", marker="o", markersize=10, alpha=0.5),
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
        Line2D([0], [0], color="k", lw=2, alpha=0.5, linestyle="--"),
        Line2D([0], [0], color="k", lw=2, linestyle=":", alpha=0.65),
    ]
    ax[1].legend(
        custom_lines,
        ["Diffuse", "Diffuse R=25", "G21 Fit", "G21 Components", "F19 R(V)=3.1"],
    )

    fig.tight_layout()

    save_fname = "diffuse_aveext"
    if args.png:
        fig.savefig(save_fname + ".png")
    elif args.pdf:
        fig.savefig(save_fname + ".pdf")
    else:
        pyplot.show()
