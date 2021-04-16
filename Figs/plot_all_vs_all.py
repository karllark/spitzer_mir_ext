# plot to plot the silicate strength versus other properties of the sightlines

import argparse
import math
import os.path

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
    alpha = np.full((n_ext), 0.0)
    alpha_unc = np.full((n_ext), 0.0)
    sil_amp = np.full((n_ext), 0.0)
    sil_amp_unc = np.full((n_ext), 0.0)
    sil_width = np.full((n_ext), 0.0)
    sil_width_unc = np.full((n_ext), 0.0)
    sil_lambda = np.full((n_ext), 0.0)
    sil_lambda_unc = np.full((n_ext), 0.0)
    sil_asym = np.full((n_ext), 0.0)
    sil_asym_unc = np.full((n_ext), 0.0)
    sil_area = np.full((n_ext), 0.0)
    sil_area_unc = np.full((n_ext), 0.0)
    sil_ceninten = np.full((n_ext), 0.0)
    sil_ceninten_unc = np.full((2, n_ext), 0.0)
    sil2_amp = np.full((n_ext), 0.0)
    sil2_amp_unc = np.full((n_ext), 0.0)
    nuv_amp = np.full((n_ext), 0.0)
    nuv_amp_unc = np.full((n_ext), 0.0)
    nuv_lambda = np.full((n_ext), 0.0)
    nuv_lambda_unc = np.full((n_ext), 0.0)
    nuv_width = np.full((n_ext), 0.0)
    nuv_width_unc = np.full((n_ext), 0.0)
    nuv_area = np.full((n_ext), 0.0)
    nuv_area_unc = np.full((n_ext), 0.0)
    nuv_ceninten = np.full((n_ext), 0.0)
    nuv_ceninten_unc = np.full((2, n_ext), 0.0)
    fuv_amp = np.full((n_ext), 0.0)
    fuv_amp_unc = np.full((n_ext), 0.0)
    avs = np.full((n_ext), 0.0)
    avs_unc = np.full((2, n_ext), 0.0)
    ebvs = np.full((n_ext), 0.0)
    ebvs_unc = np.full((n_ext), 0.0)
    rvs = np.full((n_ext), 0.0)
    rvs_unc = np.full((2, n_ext), 0.0)

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

        # print(np.corrcoef(np.stack([avs_dist.distribution, rvs_dist.distribution])))

        alpha_dist = unc.Distribution(samples[:, 1])
        silamp_dist = unc.Distribution(samples[:, 2])
        sillam_dist = unc.Distribution(samples[:, 3])
        silwid_dist = unc.Distribution(samples[:, 4])
        silasm_dist = unc.Distribution(samples[:, 5])
        silamp2_dist = unc.Distribution(samples[:, 6])

        alpha[k] = alpha_dist.pdf_mean()
        alpha_unc[k] = alpha_dist.pdf_std()

        sil_amp[k] = silamp_dist.pdf_mean()
        sil_amp_unc[k] = silamp_dist.pdf_std()
        sil_width[k] = silwid_dist.pdf_mean()
        sil_width_unc[k] = silwid_dist.pdf_std()
        sil_lambda[k] = sillam_dist.pdf_mean()
        sil_lambda_unc[k] = sillam_dist.pdf_std()
        sil_asym[k] = silasm_dist.pdf_mean()
        sil_asym_unc[k] = silasm_dist.pdf_std()

        silceninten_dist = silamp_dist / ((silwid_dist / sillam_dist) ** 2)
        sil_ci_per = silceninten_dist.pdf_percentiles([16.0, 50.0, 84.0])
        sil_ceninten[k] = sil_ci_per[1]
        sil_ceninten_unc[1, k] = sil_ci_per[2] - sil_ci_per[1]
        sil_ceninten_unc[0, k] = sil_ci_per[1] - sil_ci_per[0]

        # using C3 = a_i * lambda_i**2 * 2 to be able to use the FM area formula
        silarea_dist = math.pi * silamp_dist * (sillam_dist ** 2) / silwid_dist
        sil_area[k] = silarea_dist.pdf_mean()
        sil_area_unc[k] = silarea_dist.pdf_std()

        sil2_amp[k] = silamp2_dist.pdf_mean()
        sil2_amp_unc[k] = silamp2_dist.pdf_std()

        # get FM90 fits
        uvfname = bfile.replace(".fits", "_FM90.fits")
        if os.path.isfile(uvfname):
            cext_fm90 = ExtData(filename=uvfname)
            mcmcfile = uvfname.replace(".fits", ".h5")
            reader = emcee.backends.HDFBackend(mcmcfile)
            nsteps, nwalkers = reader.get_log_prob().shape
            samples = reader.get_chain(discard=int(mcmc_burnfrac * nsteps), flat=True)

            nuvamp_dist = unc.Distribution(samples[:, 2])
            nuvlam_dist = unc.Distribution(samples[:, 4])
            nuvwid_dist = unc.Distribution(samples[:, 5])

            fuvamp_dist = unc.Distribution(samples[:, 3])

            nuv_amp[k] = nuvamp_dist.pdf_mean()
            nuv_amp_unc[k] = nuvamp_dist.pdf_std()
            nuv_width[k] = nuvwid_dist.pdf_mean()
            nuv_width_unc[k] = nuvwid_dist.pdf_std()
            nuv_lambda[k] = nuvlam_dist.pdf_mean()
            nuv_lambda_unc[k] = nuvlam_dist.pdf_std()

            fuv_amp[k] = fuvamp_dist.pdf_mean()
            fuv_amp_unc[k] = fuvamp_dist.pdf_std()

            nuvceninten_dist = nuvamp_dist / (nuvwid_dist ** 2)
            nuv_ci_per = nuvceninten_dist.pdf_percentiles([16.0, 50.0, 84.0])
            nuv_ceninten[k] = nuv_ci_per[1]
            nuv_ceninten_unc[1, k] = nuv_ci_per[2] - nuv_ci_per[1]
            nuv_ceninten_unc[0, k] = nuv_ci_per[1] - nuv_ci_per[0]

            # using C3 = a_i * lambda_i**2 * 2 to be able to use the FM area formula
            nuvarea_dist = math.pi * nuvamp_dist / (2.0 * nuvwid_dist)
            nuv_area[k] = nuvarea_dist.pdf_mean()
            nuv_area_unc[k] = nuvarea_dist.pdf_std()

    # output some info
    a = Table()
    a["AV"] = avs
    # a["RV"] = rvs
    a["alpha"] = alpha
    a["Sil1_amp"] = sil_amp
    a["Sil1_lambda"] = sil_lambda
    a["Sil1_width"] = sil_width
    a["Sil1_asym"] = sil_asym
    # a["Sil2_amp"] = sil2_amp

    # output some info
    a_unc = Table()
    a_unc["AV"] = 0.5 * (avs_unc[0, :] + avs_unc[0, :])
    # a_unc["RV"] = 0.5 * (rvs_unc[0, :] + rvs_unc[0, :])
    a_unc["alpha"] = alpha_unc
    a_unc["Sil1_amp"] = sil_amp_unc
    a_unc["Sil1_lambda"] = sil_lambda_unc
    a_unc["Sil1_width"] = sil_width_unc
    a_unc["Sil1_asym"] = sil_asym_unc
    # a_unc["Sil2_amp"] = sil2_amp_unc

    b = Table()
    gooduv = nuv_amp > 0
    b["AV"] = avs
    b["RV"] = rvs
    b["NUV_amp"] = nuv_amp
    b["NUV_lambda"] = nuv_lambda
    b["NUV_width"] = nuv_width
    b["FUV_amp"] = fuv_amp

    b_unc = Table()
    b_unc["AV"] = 0.5 * (avs_unc[0, :] + avs_unc[0, :])
    b_unc["RV"] = 0.5 * (rvs_unc[0, :] + rvs_unc[0, :])
    b_unc["NUV_amp"] = nuv_amp_unc
    b_unc["NUV_lambda"] = nuv_lambda_unc
    b_unc["NUV_width"] = nuv_width_unc
    b_unc["FUV_amp"] = fuv_amp_unc

    # plots
    fontsize = 10

    font = {"size": fontsize}

    matplotlib.rc("font", **font)

    matplotlib.rc("lines", linewidth=1)
    matplotlib.rc("axes", linewidth=2)
    matplotlib.rc("xtick.major", width=2)
    matplotlib.rc("xtick.minor", width=2)
    matplotlib.rc("ytick.major", width=2)
    matplotlib.rc("ytick.minor", width=2)

    figsize = (11.0, 9.0)
    gkeys = a.keys()
    nkeys = len(gkeys)
    fig, ax = pyplot.subplots(nrows=nkeys, ncols=nkeys, figsize=figsize)

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

    uvdiffuse = np.logical_and(diffuse, gooduv)
    uvdense = np.logical_and(dense, gooduv)

    # plot MIR parameters
    for i, ikey in enumerate(gkeys):
        for j, jkey in enumerate(gkeys):
            if j == i:
                ax[i, j].set_axis_off()
            if j < i:
                ax[i, j].errorbar(
                    a[jkey].data[diffuse],
                    a[ikey].data[diffuse],
                    xerr=a_unc[jkey].data[diffuse],
                    yerr=a_unc[ikey].data[diffuse],
                    fmt="go",
                )
                ax[i, j].errorbar(
                    a[jkey].data[dense],
                    a[ikey].data[dense],
                    xerr=a_unc[jkey].data[dense],
                    yerr=a_unc[ikey].data[dense],
                    fmt="bo",
                    markerfacecolor="none",
                )
            if i == (nkeys - 1):
                ax[i, j].set_xlabel(jkey)
            if j == 0:
                ax[i, j].set_ylabel(ikey)

    for i, ikey in enumerate(b.keys()):
        for j, jkey in enumerate(b.keys()):
            if j > i:
                ax[i, j].errorbar(
                    b[jkey].data[uvdiffuse],
                    b[ikey].data[uvdiffuse],
                    xerr=b_unc[jkey].data[uvdiffuse],
                    yerr=b_unc[ikey].data[uvdiffuse],
                    fmt="go",
                )
                ax[i, j].errorbar(
                    b[jkey].data[uvdense],
                    b[ikey].data[uvdense],
                    xerr=b_unc[jkey].data[uvdense],
                    yerr=b_unc[ikey].data[uvdense],
                    fmt="bo",
                    markerfacecolor="none",
                )
                ax[i, j].xaxis.tick_top()
                ax[i, j].xaxis.set_label_position("top")
                ax[i, j].yaxis.tick_right()
                ax[i, j].yaxis.set_label_position("right")

            if i == 0:
                ax[i, j].set_xlabel(jkey)
            if j == len(b.keys()) - 1:
                ax[i, j].set_ylabel(ikey)

    fig.tight_layout()

    save_str = "_all_vs_all"
    if args.png:
        fig.savefig(args.filelist.replace(".dat", save_str + ".png"))
    elif args.pdf:
        fig.savefig(args.filelist.replace(".dat", save_str + ".pdf"))
    else:
        pyplot.show()
