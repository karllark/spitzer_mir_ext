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
        "XO": r"{[$\micron^{-1}$]}",
        "GAMMA": r"{[$\micron^{-1}$]}",
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

    for k, bfile in enumerate(files):
        edata = ExtData(filename=bfile.replace(".fits", "_FM90.fits"))

        spos = names[k].find("_")
        sname = names[k][:spos].upper()

        pstr = f"{sname} & "
        for k, ckey in enumerate(okeys):
            if first_line:
                hstr += fr"\colhead{phead[ckey]} & "
                hstr2 += fr"\colhead{phead2[ckey]} & "
            val, punc, munc = edata.fm90_p50_fit[ckey]
            cmval = float(mval[ckey])
            if sname == "DIFFUS":
                pstr += f"${cmval*val:.4f}^{{+{cmval*punc:.4f}}}_{{-{cmval*munc:.4f}}}$ & "
            else:
                pstr += f"${cmval*val:.3f}^{{+{cmval*punc:.3f}}}_{{-{cmval*munc:.3f}}}$ & "
        if first_line:
            first_line = False
            print(f"\\tablehead{{{hstr[:-3]}}} \\\\")
            print(f"\\{{{hstr2[:-3]}}}")
            print(f"\\startdata")
        print(f"{pstr[:-3]} \\\\")
