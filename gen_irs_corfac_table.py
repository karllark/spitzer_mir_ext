#!/usr/bin/env python
#
# Program to generate the table of IRS correction factors
#
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import argparse

import numpy as np

from getstardata import StarData

if __name__ == "__main__":

    # commandline parser
    parser = argparse.ArgumentParser()
    parser.add_argument("filelist", help="file with list of stars to use")
    args = parser.parse_args()

    filename = args.filelist

    f = open(filename, 'r')
    file_lines = list(f)
    starnames = []
    stardata = []
    bvcol = []
    for line in file_lines:
        if (line.find('#') != 0) & (len(line) > 0):
            name = line.rstrip()
            starnames.append(name)
            tstar = StarData('DAT_files/'+name+'.dat', 
                             path='/home/kgordon/Dust/Ext/')
            if 'IRS_slope' in tstar.corfac.keys():
                print('%s & %.3f & %.3f & %.3f \\\\'%(name.upper(), 
                                            tstar.corfac['IRS'], 
                                            tstar.corfac['IRS_slope'], 
                                            tstar.corfac['IRS_zerowave']))
            else:
                print('%s & %.3f & \\nodata & \\nodata \\\\'%(name.upper(), 
                                            tstar.corfac['IRS']))
