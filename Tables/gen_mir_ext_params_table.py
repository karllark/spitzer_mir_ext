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
        "SIL1_AMP": r"{$10^{-3}$ $A(\lambda)/A(V)$}",
        "SIL1_LAMBDA": r"{[$\micron$]}",
        "SIL1_WIDTH": r"{[$\micron$]}",
        "SIL2_AMP": r"{$10^{-3}$ $A(\lambda)/A(V)$}",
        "FIR_AMP": r"{$10^{-3}$ $A(\lambda)/A(V)$}",
    }
    mval = {
        "AV": 1,
        "SIL1_AMP": 1e3,
        "SIL1_LAMBDA": 1,
        "SIL1_WIDTH": 1,
        "SIL2_AMP": 1e3,
        "FIR_AMP": 1e3,
    }

    okeys = ["AV", "SIL1_AMP", "SIL1_LAMBDA", "SIL1_WIDTH", "SIL2_AMP", "FIR_AMP"]

    for line in sorted(file_lines):
        if (line.find("#") != 0) & (len(line) > 0):
            name = line.rstrip()
            edata = ExtData(filename=f"fits/{name}")

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
                else:
                    val, punc, munc = edata.p92_p50_fit[ckey]
                cmval = float(mval[ckey])
                pstr += f"${cmval*val:.3f}^{{+{cmval*punc:.3f}}}_{{-{cmval*munc:.3f}}}$ & "
            if first_line:
                first_line = False
                print(f"\\tablehead{{{hstr[:-3]}}} \\\\")
                print(f"\\{{{hstr2[:-3]}}}")
                print(f"\\startdata")
            print(f"{pstr[:-3]} \\\\")
