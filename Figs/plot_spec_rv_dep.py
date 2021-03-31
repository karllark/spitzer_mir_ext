import argparse

import emcee
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import ScalarFormatter
import astropy.units as u
from astropy import uncertainty as unc
from astropy.modeling import models, fitting

from measure_extinction.extdata import ExtData

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filelist", help="file with list of curves to plot")
    parser.add_argument(
        "--irv", help="fit versus 1/R(V) instead of R(V) - 3.1", action="store_true"
    )
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    filename = args.filelist

    # data
    f = open(filename, "r")
    file_lines = list(f)
    extfnames = []
    extnames = []
    for line in file_lines:
        if (line.find("#") != 0) & (len(line) > 0):
            name = line.rstrip()
            extfnames.append(name)
            extnames.append(name.split("_")[0])

    n_ext = len(extnames)
    avs = np.full((n_ext), 0.0)
    avs_unc = np.full((2, n_ext), 0.0)
    ebvs = np.full((n_ext), 0.0)
    ebvs_unc = np.full((n_ext), 0.0)
    rvs = np.full((n_ext), 0.0)
    rvs_unc = np.full((2, n_ext), 0.0)
    irvs = np.full((n_ext), 0.0)
    irvs_unc = np.full((2, n_ext), 0.0)

    # bands have to be treated separately as they are variable in length
    bandwaves = [1.235, 1.662, 2.159, 3.52, 4.45, 5.66, 7.67, 15.4, 23.36]
    n_band = len(bandwaves)
    bandexts = np.zeros((n_band, n_ext))
    bandexts_unc = np.zeros((n_band, n_ext))
    bandexts_npts = np.zeros((n_band, n_ext))

    for k, cname in enumerate(extfnames):

        bfile = f"fits/{cname}"
        cext = ExtData(filename=bfile)

        mcmcfile = bfile.replace(".fits", ".h5")
        reader = emcee.backends.HDFBackend(mcmcfile)
        nsteps, nwalkers = reader.get_log_prob().shape
        samples = reader.get_chain(discard=int(0.4 * nsteps), flat=True)

        avs_dist = unc.Distribution(samples[:, -1])
        av_per = avs_dist.pdf_percentiles([16.0, 50.0, 84.0])
        avs[k] = av_per[1]
        avs_unc[1, k] = av_per[2] - av_per[1]
        avs_unc[0, k] = av_per[1] - av_per[0]
        # print(avs_dist.pdf_percentiles([33., 50., 87.]))

        (indxs,) = np.where(
            (cext.waves["BAND"] > 0.4 * u.micron)
            & (cext.waves["BAND"] < 0.5 * u.micron)
        )
        ebvs_dist = unc.normal(
            cext.exts["BAND"][indxs[0]],
            std=cext.uncs["BAND"][indxs[0]],
            n_samples=avs_dist.n_samples,
        )
        ebvs[k] = ebvs_dist.pdf_mean()
        ebvs_unc[k] = ebvs_dist.pdf_std()

        rvs_dist = avs_dist / ebvs_dist
        rv_per = rvs_dist.pdf_percentiles([16.0, 50.0, 84.0])
        rvs[k] = rv_per[1]
        rvs_unc[1, k] = rv_per[2] - rv_per[1]
        rvs_unc[0, k] = rv_per[1] - rv_per[0]

        irvs_dist = ebvs_dist / avs_dist
        irv_per = irvs_dist.pdf_percentiles([16.0, 50.0, 84.0])
        irvs[k] = irv_per[1]
        irvs_unc[1, k] = irv_per[2] - irv_per[1]
        irvs_unc[0, k] = irv_per[1] - irv_per[0]

        # get the IRS data
        if k == 0:
            n_irs = len(cext.exts["IRS"])
            irswaves = cext.waves["IRS"]
            irsexts = np.zeros((n_irs, n_ext))
            irsexts_unc = np.zeros((n_irs, n_ext))
            irsexts_npts = np.zeros((n_irs, n_ext))

        irsexts[:, k] = cext.exts["IRS"]
        irsexts_unc[:, k] = cext.uncs["IRS"]
        irsexts_npts[:, k] = cext.npts["IRS"]

        # get the A(lambda)/A(V) at certain wavelengths
        cwaves = cext.waves["BAND"]
        for j, wave in enumerate(bandwaves):
            indx = np.abs(cwaves.value - wave).argmin()
            if (np.abs(cwaves[indx].value - wave)) < 0.01:
                bandexts[j, k] = cext.exts["BAND"][indx]
                bandexts_unc[j, k] = cext.uncs["BAND"][indx]
                bandexts_npts[j, k] = cext.npts["BAND"][indx]

    # for every wavelength, fit a straight line through the A(lambda)/A(V) vs. 1/R(V) data
    fit = fitting.LinearLSQFitter()
    line_func = models.Linear1D()
    slopes = []
    intercepts = []
    stds = []
    wave_list = irswaves
    for j, wave in enumerate(wave_list):
        mask = irsexts_npts[j] > 0
        if args.irv:
            xvals = irvs[mask]
        else:
            xvals = rvs[mask] - 3.1
        yvals = (irsexts[j][mask] / avs[mask]) + 1.0
        yuncs = irsexts_unc[j][mask] / avs[mask]
        fitted_line = fit(line_func, xvals, yvals, weights=1.0 / yuncs)
        std = np.sqrt(np.sum((fitted_line(xvals) - yvals) ** 2) / len(xvals))
        slopes.append(fitted_line.slope.value)
        intercepts.append(fitted_line.intercept.value)
        stds.append(std)

    bslopes = []
    bintercepts = []
    bstds = []
    for j, wave in enumerate(bandwaves):
        mask = bandexts_npts[j] > 0
        if args.irv:
            xvals = irvs[mask]
        else:
            xvals = rvs[mask] - 3.1
        yvals = (bandexts[j][mask] / avs[mask]) + 1.0
        yuncs = bandexts_unc[j][mask] / avs[mask]
        fitted_line = fit(line_func, xvals, yvals, weights=1.0 / yuncs)
        std = np.sqrt(np.sum((fitted_line(xvals) - yvals) ** 2) / len(xvals))
        bslopes.append(fitted_line.slope.value)
        bintercepts.append(fitted_line.intercept.value)
        bstds.append(std)

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

    # plot the slopes and intercepts vs. wavelength
    figsize = (10.0, 8.0)
    fig, ax = plt.subplots(3, figsize=figsize, sharex=True)
    ax[0].plot(irswaves, slopes, "bo", markersize=2)
    ax[0].plot(bandwaves, bslopes, "go")
    ax[1].plot(irswaves, intercepts, "bo", markersize=2)
    ax[1].plot(bandwaves, bintercepts, "go")
    ax[2].plot(irswaves, stds, "bo", markersize=2)
    ax[2].plot(bandwaves, bstds, "go")

    # add CCM89 RV dependent relation
    if args.irv:
        waves = [0.7, 0.9, 1.25, 1.6, 2.2, 3.4]
        intercepts = [0.8686, 0.68, 0.4008, 0.2693, 0.1615, 0.08]
        slopes = [-0.366, -0.6239, -0.3679, -0.2473, -0.1483, -0.0734]
        ax[0].scatter(waves, slopes, s=5)
        ax[1].scatter(waves, intercepts, s=5)

    plt.xlabel(r"$\lambda$ [$\mu m$]")
    ax[0].set_ylabel("slopes")
    ax[1].set_ylabel("intercepts")
    ax[2].set_ylabel("stddevs")
    # plt.subplots_adjust(hspace=0)
    # plt.savefig(outpath + "RV_slope_inter.pdf", bbox_inches="tight")
    ax[0].set_xlim(1.0, 35.0)
    ax[0].set_xscale("log")
    if args.irv:
        ax[1].set_ylim(-0.3, 0.3)
        ax[0].set_ylim(-0.5, 0.5)
    else:
        ax[1].set_ylim(0.0, 0.3)
        ax[0].set_ylim(-0.05, 0.05)

    ax[0].plot([0.5, 40.0], [0.0, 0.0], "k:", linewidth=3, alpha=0.5)
    ax[1].plot([0.5, 40.0], [0.0, 0.0], "k:", linewidth=3, alpha=0.5)

    ax[0].tick_params("both", length=10, width=2, which="major")
    ax[0].tick_params("both", length=5, width=1, which="minor")
    ax[1].tick_params("both", length=10, width=2, which="major")
    ax[1].tick_params("both", length=5, width=1, which="minor")

    ax[1].xaxis.set_major_formatter(ScalarFormatter())
    ax[1].xaxis.set_minor_formatter(ScalarFormatter())

    fig.tight_layout()

    save_str = "_spec_rvdep"
    if args.irv:
        save_str += "_irvs"
    if args.png:
        fig.savefig(args.filelist.replace(".dat", save_str + ".png"))
    elif args.pdf:
        fig.savefig(args.filelist.replace(".dat", save_str + ".pdf"))
    else:
        plt.show()
