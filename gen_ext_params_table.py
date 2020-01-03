#!/usr/bin/env python
#
# Program to generate the table of extinction parameters
#
import argparse

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

    phead = {
        "AV": "{A(V)}",
        "BKG_amp": "{BKG A}",
        "FUV_amp": "{FUV A}",
        "NUV_amp": "{NUV A}",
        "NUV_lambda": r"{NUV $\lambda_o$}",
        "NUV_width": r"{NUV $\gamma$}",
        "SIL1_amp": "{SIL1 A}",
        "SIL1_lambda": r"{SIL1 $\lambda_o$}",
        "SIL1_width": r"{SIL1 $\gamma$}",
        "SIL2_amp": "{SIL2 A}",
        "FIR_amp": "{FIR A}",
    }

    for line in file_lines:
        if (line.find("#") != 0) & (len(line) > 0):
            name = line.rstrip()
            edata = ExtData(filename=f"fits/{name}")

            spos = name.find("_")
            sname = name[:spos]

            pstr = f"{sname} & "
            for ckey in edata.p92_p50_fit.keys():
                if first_line:
                    hstr += fr"\colhead{phead[ckey]} & "
                val, punc, munc = edata.p92_p50_fit[ckey]
                pstr += f"${val:.3f}^{{+{punc:.3f}}}_{{-{munc:.3f}}}$ & "
            if first_line:
                first_line = False
                print(f"\\tablehead{{{hstr[:-3]}}}")
                print(f"\\startdata")
            print(f"{pstr[:-3]} \\\\")
