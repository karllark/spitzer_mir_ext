import argparse

import emcee
import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib
import astropy.units as u
from astropy import uncertainty as unc
from astropy.modeling import models, fitting

from measure_extinction.extdata import ExtData

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filelist", help="file with list of curves to plot")
    parser.add_argument("--band", help="band for plot", default="U")
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
    exts = np.zeros(n_ext)
    exts_unc = np.zeros(n_ext)

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

        (indxs,) = np.where(cext.names["BAND"] == args.band)
        exts[k] = (cext.exts["BAND"][indxs[0]] / avs[k]) + 1.0
        exts_unc[k] = cext.uncs["BAND"][indxs[0]] / avs[k]

    # fit
    fit = fitting.LinearLSQFitter()
    line_func = models.Linear1D()
    fitted_line = fit(line_func, (rvs - 3.1), exts)
    print(fitted_line)

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

    figsize = (5.5, 5.0)
    fig, ax = pyplot.subplots(figsize=figsize)

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

    # R(V) versus A(V)
    ax.errorbar(
        rvs[diffuse] - 3.1,
        exts[diffuse],
        xerr=rvs_unc[:, diffuse],
        yerr=exts_unc[diffuse],
        fmt="go",
        label="diffuse",
    )
    ax.errorbar(
        rvs[dense] - 3.1,
        exts[dense],
        xerr=irvs_unc[:, dense],
        yerr=exts_unc[dense],
        fmt="bo",
        markerfacecolor="none",
        label="dense",
    )

    mrvs = np.array([-0.75, 1.5])
    ax.plot(mrvs, fitted_line(mrvs), "k:")

    ax.set_xlabel(r"$R(V) - 3.1$")
    ax.set_ylabel(rf"$A({args.band})/A(V)$")
    ax.tick_params("both", length=10, width=2, which="major")
    ax.tick_params("both", length=5, width=1, which="minor")

    ax.legend()

    fig.tight_layout()

    save_str = f"_rvdep{args.band}"
    if args.png:
        fig.savefig(args.filelist.replace(".dat", save_str + ".png"))
    elif args.pdf:
        fig.savefig(args.filelist.replace(".dat", save_str + ".pdf"))
    else:
        pyplot.show()
