#!/usr/bin/env python
#
# Program to generate the table of extinction parameters
#
import argparse

import emcee
# from astropy import uncertainty as unc

from measure_extinction.extdata import ExtData

if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("filelist", help="file with list of stars to use")
    args = parser.parse_args()

    filename = args.filelist
    f = open(filename, "r")
    file_lines = list(f)

    files = []
    names = []
    for line in sorted(file_lines):
        if (line.find("#") != 0) & (len(line) > 0):
            name = line.rstrip()
            names.append(name)
            bfile = f"fits/{name}"
            files.append(bfile)
    # add diffuse average
    names.append("diffuse")
    files.append("data/all_ext_18feb20_diffuse_ave_POWLAW2DRUDE.fits")

    first_line = True
    hstr = r"\colhead{name} & "
    hstr2 = r" & "

    phead = {
        "AV": "{A(V)}",
        "SIL1_AMP": "{SIL1 A}",
        "SIL1_CENTER": r"{SIL1 $\lambda_o$}",
        "SIL1_FWHM": r"{SIL1 $\gamma$}",
        "SIL1_ASYM": "{SIL1 e}",
        "SIL2_AMP": "{SIL2 A}",
        "SIL2_CENTER": r"{SIL2 $\lambda_o$}",
        "SIL2_FWHM": r"{SIL2 $\gamma$}",
        "SIL2_ASYM": "{SIL2 e}",
    }
    phead2 = {
        "AV": "{[mag]}",
        "SIL1_AMP": r"{$10^{-2}$ $A(\lambda)/A(V)$}",
        "SIL1_CENTER": r"{[$\micron$]}",
        "SIL1_FWHM": r"{[$\micron$]}",
        "SIL1_ASYM": r"{}",
        "SIL2_AMP": r"{$10^{-2}$ $A(\lambda)/A(V)$}",
        "SIL2_CENTER": r"{[$\micron$]}",
        "SIL2_FWHM": r"{[$\micron$]}",
        "SIL2_ASYM": r"{}",
    }
    mval = {
        "AV": 1,
        "SIL1_AMP": 1e2,
        "SIL1_CENTER": 1,
        "SIL1_FWHM": 1,
        "SIL1_ASYM": 1,
        "SIL2_AMP": 1e2,
        "SIL2_CENTER": 1,
        "SIL2_FWHM": 1,
        "SIL2_ASYM": 1,
    }

    # fmt: off
    okeys = ["AV", "SIL1_AMP", "SIL1_CENTER", "SIL1_FWHM", "SIL1_ASYM",
             "SIL2_AMP", "SIL2_CENTER", "SIL2_FWHM", "SIL2_ASYM"]
    # fmt: on

    mcmc_burnfrac = 0.4
    for k, bfile in enumerate(files):
        edata = ExtData(filename=bfile)

        mcmcfile = bfile.replace(".fits", ".h5")
        reader = emcee.backends.HDFBackend(mcmcfile)
        nsteps, nwalkers = reader.get_log_prob().shape
        samples = reader.get_chain(discard=int(mcmc_burnfrac * nsteps), flat=True)

        spos = names[k].find("_")
        sname = names[k][:spos].upper()

        pstr = f"{sname} & "
        for k, ckey in enumerate(okeys):
            if first_line:
                hstr += fr"\colhead{phead[ckey]} & "
                hstr2 += fr"\colhead{phead2[ckey]} & "
                # hstr2 += fr"\colhead{{{mval[k]:.1f}}} & "
            if ckey == "AV":
                if hasattr(edata, "columns_p50_fit"):
                    val, punc, munc = edata.columns_p50_fit[ckey]
                else:
                    val, punc, munc = (1.0, 0.0, 0.0)
            else:
                val, punc, munc = edata.g21_p50_fit[ckey]
            cmval = float(mval[ckey])
            pstr += (
                f"${cmval*val:.2f}^{{+{cmval*punc:.2f}}}_{{-{cmval*munc:.2f}}}$ & "
            )
        if first_line:
            first_line = False
            print(f"\\tablehead{{{hstr[:-3]}}} \\\\")
            print(f"\\{{{hstr2[:-3]}}}")
            print(f"\\startdata")
        print(f"{pstr[:-3]} \\\\")
