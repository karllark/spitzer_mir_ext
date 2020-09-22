import argparse

import emcee
import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib
import astropy.units as u
from astropy.table import Table
from astropy import uncertainty as unc

from measure_extinction.extdata import ExtData

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filelist", help="file with list of curves to plot")
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    filename = args.filelist

    mcmc_burnfrac = 0.4

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

    # Amplitude of the feature when x = x_o (FM90) or lambda = lambda_o (P92)
    # FM90: I(x_o) = C3/(2*gamma**2)
    # P92: I(lambda_o) = a_i/(2 + b_i) = a_i/(gamma_i**2/lambda_i**2)
    # C3 = a_i * lambda_i**2 * 2.

    n_ext = len(extnames)
    avs = np.full((n_ext), 0.0)
    avs_unc = np.full((2, n_ext), 0.0)
    ebvs = np.full((n_ext), 0.0)
    ebvs_unc = np.full((n_ext), 0.0)
    ekvs = np.full((n_ext), 0.0)
    ekvs_unc = np.full((n_ext), 0.0)
    eivs = np.full((n_ext), 0.0)
    eivs_unc = np.full((n_ext), 0.0)
    rvs = np.full((n_ext), 0.0)
    rvs_unc = np.full((2, n_ext), 0.0)
    rks = np.full((n_ext), 0.0)
    rks_unc = np.full((2, n_ext), 0.0)
    ris = np.full((n_ext), 0.0)
    ris_unc = np.full((2, n_ext), 0.0)

    for k, cname in enumerate(extfnames):

        # get P92 fits
        bfile = f"fits/{cname}"
        cext = ExtData(filename=bfile)

        mcmcfile = bfile.replace(".fits", ".h5")
        reader = emcee.backends.HDFBackend(mcmcfile)
        nsteps, nwalkers = reader.get_log_prob().shape
        samples = reader.get_chain(discard=int(mcmc_burnfrac * nsteps), flat=True)

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

        (indxs,) = np.where(
            (cext.waves["BAND"] > 2.1 * u.micron)
            & (cext.waves["BAND"] < 2.3 * u.micron)
        )
        # print(cext.waves["BAND"][indxs[0]])
        ekvs_dist = unc.normal(
            cext.exts["BAND"][indxs[0]],
            std=cext.uncs["BAND"][indxs[0]],
            n_samples=avs_dist.n_samples,
        )
        ekvs[k] = ekvs_dist.pdf_mean()
        ekvs_unc[k] = ekvs_dist.pdf_std()

        rks_dist = avs_dist / ekvs_dist
        rk_per = rks_dist.pdf_percentiles([16.0, 50.0, 84.0])
        rks[k] = rk_per[1]
        rks_unc[1, k] = rk_per[2] - rk_per[1]
        rks_unc[0, k] = rk_per[1] - rk_per[0]

        (indxs,) = np.where(
            (cext.waves["BAND"] > 5.0 * u.micron)
            & (cext.waves["BAND"] < 6.0 * u.micron)
        )
        # print(cext.waves["BAND"][indxs[0]])
        eivs_dist = unc.normal(
            cext.exts["BAND"][indxs[0]],
            std=cext.uncs["BAND"][indxs[0]],
            n_samples=avs_dist.n_samples,
        )
        eivs[k] = ekvs_dist.pdf_mean()
        eivs_unc[k] = ekvs_dist.pdf_std()

        ris_dist = avs_dist / eivs_dist
        ri_per = ris_dist.pdf_percentiles([16.0, 50.0, 84.0])
        ris[k] = ri_per[1]
        ris_unc[1, k] = ri_per[2] - ri_per[1]
        ris_unc[0, k] = ri_per[1] - ri_per[0]

    # output some info
    a = Table()
    a["name"] = extnames
    a["AV"] = avs
    a["RV"] = rvs
    # print(a)

    # compute mean rk
    gvals = rks > -1.15
    mean_rk = np.average(rks[gvals], weights=0.5 * np.sum(rks_unc[:, gvals], axis=0))

    gvals = ris > -2.0
    mean_ri = np.average(ris[gvals], weights=0.5 * np.sum(ris_unc[:, gvals], axis=0))
    print(mean_rk, mean_ri)

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

    figsize = (11.0, 5.0)
    fig, fax = pyplot.subplots(ncols=2, figsize=figsize)

    diffuse = []
    for tname in extnames:
        if tname == "hd283809":
            diffuse.append(False)
        elif tname == "hd029647":
            diffuse.append(False)
        else:
            diffuse.append(True)
    diffuse = np.array(diffuse)
    dense = ~diffuse

    # A(V)/E(K-V) versus A(V)
    ax = fax[0]
    ax.errorbar(
        avs[diffuse],
        rks[diffuse],
        xerr=avs_unc[:, diffuse],
        yerr=rks_unc[:, diffuse],
        fmt="go",
        label="diffuse",
    )
    ax.errorbar(
        avs[dense],
        rks[dense],
        xerr=avs_unc[:, dense],
        yerr=rks_unc[:, dense],
        fmt="bo",
        markerfacecolor="none",
        label="dense",
    )
    ax.axhline(mean_rk, ls=":", lw=2, alpha=0.5, color="k", label="clipped average")
    ax.set_ylabel(r"$A(V)/E(K-V)$")
    ax.set_xlabel(r"$A(V)$")
    ax.tick_params("both", length=10, width=2, which="major")
    ax.tick_params("both", length=5, width=1, which="minor")
    ax.legend()

    # A(V)/E(K-V) versus R(V)
    ax = fax[1]
    ax.errorbar(
        avs[diffuse],
        ris[diffuse],
        xerr=avs_unc[:, diffuse],
        yerr=ris_unc[:, diffuse],
        fmt="go",
        label="diffuse",
    )
    ax.errorbar(
        avs[dense],
        ris[dense],
        xerr=avs_unc[:, dense],
        yerr=ris_unc[:, dense],
        fmt="bo",
        markerfacecolor="none",
        label="dense",
    )
    ax.axhline(mean_ri, ls=":", lw=2, alpha=0.5, color="k", label="average")
    ax.set_ylabel(r"$A(V)/E(I3-V)$")
    ax.set_xlabel(r"$A(V)$")
    ax.tick_params("both", length=10, width=2, which="major")
    ax.tick_params("both", length=5, width=1, which="minor")
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.legend(loc="lower right")

    fig.tight_layout()

    save_str = "_avekv"
    if args.png:
        fig.savefig(args.filelist.replace(".dat", save_str + ".png"))
    elif args.pdf:
        fig.savefig(args.filelist.replace(".dat", save_str + ".pdf"))
    else:
        pyplot.show()
