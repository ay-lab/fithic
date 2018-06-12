#!/usr/bin/env python

"""
Created: February 06 2018

Code for creating fragments to be used in FitHiC

Author: Sourya Bhattacharyya
Vijay-Ay lab, LJI
"""
## Modified by Arya Kaul 2018-05-30 to use gzip, python3, and argparse
## Modified by Arya Kaul 2018-05-31 to fix last bin bug 


import sys
import os
import argparse 
import gzip
#===============================================
def parseargs(arguments):
    parser = argparse.ArgumentParser(description="Check help flag")
    parser.add_argument("--chrLens", help="Chromosome lengths file. In format, 'chrNUM\tlength'", required=True)
    parser.add_argument("--outFile", help="Output file for storing fragment file", required=True)
    parser.add_argument("--resolution", help="Resolution of dataset being analyzed", type=int, required=True)
    return parser.parse_args()


def main():
    args = parseargs(sys.argv[1:])
    global ChrSizeFile
    global OutFile
    global BinSize
    ChrSizeFile = args.chrLens
    OutFile = args.outFile
    BinSize = args.resolution

    if ChrSizeFile is not None:
        ChrSizeFile = os.path.realpath(ChrSizeFile)
    else:
        sys.exit("Input reference chromosome size file is not specified - quit !!")

    print ('ChrSizeFile: %s' % ChrSizeFile)

    if OutFile is not None:
        OutFile = os.path.realpath(OutFile)
    else:
        sys.exit("Output file for storing the fragments is not specified - quit !!") 

    print('OutFile: %s'% OutFile)
    print('BinSize: %s' % BinSize)

    # open the output fragment file
    fp_out = gzip.open(OutFile, 'wt')

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
                #interval_end = (int(curr_chr_size / BinSize)) * BinSize
                interval_end = (int((curr_chr_size+BinSize) / BinSize)) * BinSize
            for val in range(0, interval_end, BinSize):
                curr_start = val
                curr_mid = val + BinSize / 2
                if (chr_count > 1) or (val > 0):
                    fp_out.write('\n')
                fp_out.write(str(curr_chr) + '\t' + str(curr_start) + '\t' + str(int(curr_mid)) + '\t' + str(1) + '\t' + str(1))    

    # close the output fragment file
    fp_out.close()

#===============================================
if __name__ == "__main__":
    main()

