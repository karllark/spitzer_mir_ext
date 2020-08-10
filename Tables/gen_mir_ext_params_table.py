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
    parser.add_argument(
        "--sil", action="store_true", help="file with list of stars to use"
    )
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

    if args.sil:
        phead = {
            "SIL1_AMP": r"{$S_1 \times 100$}",
            "SIL1_CENTER": r"{$\lambda_{o1}$}",
            "SIL1_FWHM": r"{$\gamma_{o1}$}",
            "SIL1_ASYM": r"{$a_1$}",
            "SIL2_AMP": r"{$S_2 \times 100$}",
            "SIL2_CENTER": r"{$\lambda_{o2}$}",
            "SIL2_FWHM": r"{$\gamma_{o2}$}",
            "SIL2_ASYM": r"{$a_2$}",
        }
        phead2 = {
            "SIL1_AMP": r"{$A(\lambda)/A(V)$}",
            "SIL1_CENTER": r"{[$\micron$]}",
            "SIL1_FWHM": r"{[$\micron$]}",
            "SIL1_ASYM": r"{}",
            "SIL2_AMP": r"{$A(\lambda)/A(V)$}",
            "SIL2_CENTER": r"{[$\micron$]}",
            "SIL2_FWHM": r"{[$\micron$]}",
            "SIL2_ASYM": r"{}",
        }
        mval = {
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
        okeys = ["SIL1_AMP", "SIL1_CENTER", "SIL1_FWHM", "SIL1_ASYM",
                 "SIL2_AMP", "SIL2_CENTER", "SIL2_FWHM", "SIL2_ASYM"]
        # fmt: on

    else:

        phead = {
            "AV": r"{$A(V)$}",
            "SCALE": r"{$B$}",
            "ALPHA": r"{$\alpha$}"
        }
        phead2 = {
            "AV": "{[mag]}",
            "SCALE": r"{$A(\lambda)/A(V)$}",
            "ALPHA": {}
        }
        mval = {
            "AV": 1,
            "SCALE": 1,
            "ALPHA": 1
        }

        okeys = ["AV", "SCALE", "ALPHA"]

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
            if (punc == 0.0) & (munc == 0.0):
                pstr += (
                    f"${cmval*val:.2f}$ & "
                )
            else:
                if sname == "DIFFUS":
                    pstr += (
                        f"${cmval*val:.3f}^{{+{cmval*punc:.3f}}}_{{-{cmval*munc:.3f}}}$ & "
                    )
                else:
                    pstr += (
                        f"${cmval*val:.2f}^{{+{cmval*punc:.2f}}}_{{-{cmval*munc:.2f}}}$ & "
                    )
        if first_line:
            first_line = False
            print(f"\\tablehead{{{hstr[:-3]}}} \\\\")
            print(f"\\{{{hstr2[:-3]}}}")
            print(f"\\startdata")
        print(f"{pstr[:-3]} \\\\")
