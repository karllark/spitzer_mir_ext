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
        "C1": r"{$C1$}",
        "C2": r"{$C2$}",
        "C3": r"{$C3$}",
        "C4": r"{$C4$}",
        "XO": r"{$x_o$}",
        "GAMMA": r"{$\gamma$}",
    }
    phead2 = {
        "C1": r"{[$A(\lambda)/A(V)$]}",
        "C2": r"{[$A(\lambda)/A(V)$]}",
        "C3": r"{[$A(\lambda)/A(V)$]}",
        "C4": r"{[$A(\lambda)/A(V)$]}",
        "XO": r"{[$\micron$]}",
        "GAMMA": r"{[$\micron$]}",
    }
    mval = {
        "C1": 1,
        "C2": 1,
        "C3": 1,
        "C4": 1,
        "XO": 1,
        "GAMMA": 1,
    }

    okeys = ["C1", "C2", "C3", "C4", "XO", "GAMMA"]

    for line in sorted(file_lines):
        if (line.find("#") != 0) & (len(line) > 0):
            name = line.rstrip()
            bfilename = f"fits_good/{name}"
            edata = ExtData(filename=bfilename.replace(".fits", "_fm90.fits"))

            spos = name.find("_")
            sname = name[:spos].upper()

            pstr = f"{sname} & "
            for k, ckey in enumerate(okeys):
                if first_line:
                    hstr += fr"\colhead{phead[ckey]} & "
                    hstr2 += fr"\colhead{phead2[ckey]} & "
                val, punc, munc = edata.fm90_p50_fit[ckey]
                cmval = float(mval[ckey])
                pstr += f"${cmval*val:.3f}^{{+{cmval*punc:.3f}}}_{{-{cmval*munc:.3f}}}$ & "
            if first_line:
                first_line = False
                print(f"\\tablehead{{{hstr[:-3]}}} \\\\")
                print(f"\\{{{hstr2[:-3]}}}")
                print(f"\\startdata")
            print(f"{pstr[:-3]} \\\\")
