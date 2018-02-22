#!/usr/bin/env python

"""
Created: February 06 2018

Code for creating fragments to be used in FitHiC

Author: Sourya Bhattacharyya
Vijay-Ay lab, LJI
"""

import sys
import os
from optparse import OptionParser

#===============================================
def main():
    parser = OptionParser()

    # the first option parses input chromosome size file
    parser.add_option("--ChrSizeFile", dest="ChrSizeFile", help="Input file containing size of the reference chromosome.")

    # second option takes the output file storing fragment information
    parser.add_option("--OutFile", dest="OutFile", help="Output file for storing the fragments.")

    # third option is the bin size 
    parser.add_option("--binsize", dest="binsize", type="int", help="Bin size employed for fragment file. DEFAULT 10000 (10 Kb resolution).")

    parser.set_defaults(ChrSizeFile=None, OutFile=None, binsize=10000)
    (options, args) = parser.parse_args()

    global ChrSizeFile
    global OutFile
    global BinSize

    if options.ChrSizeFile is not None:
        ChrSizeFile = os.path.realpath(options.ChrSizeFile)
    else:
        sys.exit("Input reference chromosome size file is not specified - quit !!")

    print 'ChrSizeFile: ', ChrSizeFile

    if options.OutFile is not None:
        OutFile = os.path.realpath(options.OutFile)
    else:
        sys.exit("Output file for storing the fragments is not specified - quit !!") 

    print 'OutFile: ', OutFile

    BinSize = int(options.binsize)

    print 'BinSize: ', BinSize

    # open the output fragment file
    fp_out = open(OutFile, 'w')

    # read the file containing chromosome size
    chr_count = 0
    with open(ChrSizeFile, 'r') as fp:
        for line in fp:
            if (line == ""):
                continue
            # individual line contains one chromosome and corresponding size
            linecontents = (line.rstrip()).split()
            curr_chr = linecontents[0]
            curr_chr_size = int(linecontents[1])
            chr_count = chr_count + 1
            if ((curr_chr_size % BinSize) == 0):
                interval_end = curr_chr_size
            else:
                interval_end = (int(curr_chr_size / BinSize)) * BinSize
            for val in range(0, interval_end, BinSize):
                curr_start = val
                curr_mid = val + BinSize / 2
                if (chr_count > 1) or (val > 0):
                    fp_out.write('\n')
                fp_out.write(str(curr_chr) + '\t' + str(curr_start) + '\t' + str(curr_mid) + '\t' + str(1) + '\t' + str(1))    

    # close the output fragment file
    fp_out.close()

#===============================================
if __name__ == "__main__":
    main()

