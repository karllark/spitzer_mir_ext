#!/usr/bin/env python
#
# Program to generate the table of extinction parameters
#
import argparse
import numpy as np

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
    mval = [1.0, 1.0, 1.0, 100.0, 10.0, 10.0, 1000.0, 1.0, 1.0, 1000.0, 1000.0]
    mval = np.full((len(phead)), 1.0)

    okeys = ["AV", "SIL1_AMP", "SIL1_LAMBDA", "SIL1_WIDTH",
             "SIL2_AMP", "FIR_AMP"]

    for line in file_lines:
        if (line.find("#") != 0) & (len(line) > 0):
            name = line.rstrip()
            edata = ExtData(filename=f"fits/{name}")

            spos = name.find("_")
            sname = name[:spos]

            pstr = f"{sname} & "
            for k, ckey in enumerate(okeys):
                if first_line:
                    hstr += fr"\colhead{phead[ckey]} & "
                    hstr2 += fr"\colhead{{{mval[k]:.1f}}} & "
                if ckey == "AV":
                    val, punc, munc = edata.columns_p50_fit[ckey]
                else:
                    val, punc, munc = edata.p92_p50_fit[ckey]
                pstr += f"${mval[k]*val:.3f}^{{+{mval[k]*punc:.3f}}}_{{-{mval[k]*munc:.3f}}}$ & "
            if first_line:
                first_line = False
                print(f"\\tablehead{{{hstr[:-3]}}}")
                print(f"\\tablehead{{{hstr2[:-3]}}}")
                print(f"\\startdata")
            print(f"{pstr[:-3]} \\\\")
