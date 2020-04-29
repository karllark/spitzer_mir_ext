#!/usr/bin/env python
#
# Program to generate the table of extinction parameters
#
import argparse

import emcee
from astropy import uncertainty as unc

from measure_extinction.extdata import ExtData

if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("filelist", help="file with list of stars to use")
    args = parser.parse_args()

    filename = args.filelist
    f = open(filename, "r")
    file_lines = list(f)

    first_line = True
    hstr = r"\colhead{name} & "
    hstr2 = r" & "

    phead = {
        "AV": "{A(V)}",
        "BKG_AMP": "{BKG A}",
        "FUV_AMP": "{FUV A}",
        "NUV_AMP": "{NUV A}",
        "NUV_LAMBDA": r"{NUV $\lambda_o$}",
        "NUV_WIDTH": r"{NUV $\gamma$}",
        "SIL1_AMP": "{SIL1 A}",
        "SIL1_LAMBDA": r"{SIL1 $\lambda_o$}",
        "SIL1_WIDTH": r"{SIL1 $\gamma$}",
        "SIL2_AMP": "{SIL2 A}",
        "FIR_AMP": "{FIR A}",
    }
    phead2 = {
        "AV": "{[mag]}",
        "SIL1_AMP": r"{$10^{-2}$ $A(\lambda)/A(V)$}",
        "SIL1_LAMBDA": r"{[$\micron$]}",
        "SIL1_WIDTH": r"{[$\micron$]}",
        "SIL2_AMP": r"{$10^{-3}$ $A(\lambda)/A(V)$}",
        "FIR_AMP": r"{$10^{-3}$ $A(\lambda)/A(V)$}",
    }
    mval = {
        "AV": 1,
        "SIL1_AMP": 1e2,
        "SIL1_LAMBDA": 1,
        "SIL1_WIDTH": 1,
        "SIL2_AMP": 1e3,
        "FIR_AMP": 1e3,
    }

    okeys = ["AV", "SIL1_AMP", "SIL1_LAMBDA", "SIL1_WIDTH", "SIL2_AMP", "FIR_AMP"]

    mcmc_burnfrac = 0.4
    for line in sorted(file_lines):
        if (line.find("#") != 0) & (len(line) > 0):
            name = line.rstrip()
            bfile = f"fits/{name}"
            edata = ExtData(filename=bfile)

            mcmcfile = bfile.replace(".fits", ".h5")
            reader = emcee.backends.HDFBackend(mcmcfile)
            nsteps, nwalkers = reader.get_log_prob().shape
            samples = reader.get_chain(discard=int(mcmc_burnfrac * nsteps), flat=True)

            silamp_dist = unc.Distribution(samples[:, 5])
            sillam_dist = unc.Distribution(samples[:, 6])
            silwid_dist = unc.Distribution(samples[:, 7])
            silceninten_dist = silamp_dist / ((silwid_dist / sillam_dist) ** 2)
            sil_ci_per = silceninten_dist.pdf_percentiles([16.0, 50.0, 84.0])
            sil_ci_vals = (
                sil_ci_per[1],
                sil_ci_per[2] - sil_ci_per[1],
                sil_ci_per[1] - sil_ci_per[0],
            )

            spos = name.find("_")
            sname = name[:spos].upper()

            pstr = f"{sname} & "
            for k, ckey in enumerate(okeys):
                if first_line:
                    hstr += fr"\colhead{phead[ckey]} & "
                    hstr2 += fr"\colhead{phead2[ckey]} & "
                    # hstr2 += fr"\colhead{{{mval[k]:.1f}}} & "
                if ckey == "AV":
                    val, punc, munc = edata.columns_p50_fit[ckey]
                elif ckey == "SIL1_AMP":
                    val, punc, munc = sil_ci_vals
                else:
                    val, punc, munc = edata.p92_p50_fit[ckey]
                cmval = float(mval[ckey])
                pstr += (
                    f"${cmval*val:.3f}^{{+{cmval*punc:.3f}}}_{{-{cmval*munc:.3f}}}$ & "
                )
            if first_line:
                first_line = False
                print(f"\\tablehead{{{hstr[:-3]}}} \\\\")
                print(f"\\{{{hstr2[:-3]}}}")
                print(f"\\startdata")
            print(f"{pstr[:-3]} \\\\")
