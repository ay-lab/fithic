#!/usr/bin/env python3

"""
Created: August 15, 2022

New Code for converting Juicer generated HiC files to FitHiC input format

Author: Sourya Bhattacharyya
Vijay-Ay lab, LJI
"""

import numpy as np
import hicstraw
import argparse 
import sys
import os

#===============================================
def parseargs(arguments):
    parser = argparse.ArgumentParser(description="Check help flag")
    parser.add_argument("--HiCFile", help="Input .hic file. Mandatory parameter.", required=True)
    parser.add_argument("--CHR1", help="Chromosome 1. Mandatory parameter.", required=True)
    parser.add_argument("--CHR2", help="Chromosome 2. Mandatory parameter.", required=True)
    parser.add_argument("--resolution", help="Resolution of contact matrix (in bp) to be dumped. Mandatory parameter.", type=int, required=True)
    parser.add_argument("--datatype", help="Type of contact to be dumped. Available options are 'observed' (default) or 'oe' (observed/expected).", default='observed')
    parser.add_argument("--Norm", help="Type of normalization. Available options: 'NONE', 'VC', 'VC_SQRT', 'KR' (default), 'SCALE'.", default='KR')
    parser.add_argument("--outFile", help="Output file for storing contact counts", required=True)    
    return parser.parse_args()


def main():
    args = parseargs(sys.argv[1:])

    global HiCFile
    global resolution
    global OutFile
    global datatype

    HiCFile = args.HiCFile
    if HiCFile is not None:
        HiCFile = os.path.realpath(HiCFile)
    else:
        sys.exit("Input Hi-C file is not specified - quit !!")
    print ('HiCFile: %s' % HiCFile)

    CHR1 = args.CHR1
    if CHR1 is None:        
        sys.exit("CHR1 is not specified - quit !!")
    print ('CHR1: %s' % CHR1)
    chr1name="chr" + str(CHR1)

    CHR2 = args.CHR2
    if CHR2 is None:        
        sys.exit("CHR2 is not specified - quit !!")
    print ('CHR2: %s' % CHR2)
    chr2name="chr" + str(CHR2)
    
    if args.resolution is None:
        sys.exit("Input Hi-C file is not specified - quit !!")
    else:
        resolution = int(args.resolution)
        print ('resolution: %s' % resolution)

    OutFile = args.outFile
    if OutFile is not None:
        OutFile = os.path.realpath(OutFile)
    else:
        sys.exit("Output file for storing the contacts is not specified - quit !!") 
    print('OutFile: %s'% OutFile)

    datatype = args.datatype
    print('datatype: %s'% datatype)

    Norm = args.Norm
    print('Norm: %s'% Norm)

    hic = hicstraw.HiCFile(HiCFile)
    chrom_list = hic.getChromosomes()
    print("list of chromosomes within the HiC file: " + ' '.join(str(x) for x in chrom_list))

    genomeid = hic.getGenomeID()
    print("list of genome ID within the HiC file: " + ' '.join(str(x) for x in genomeid))

    reslist = hic.getResolutions()
    print("list of resolutions within the HiC file: " + ' '.join(str(x) for x in reslist))

    if resolution not in reslist:
        sys.exit("Input resolution is not included in the given .hic file - quit !!") 

    with open(OutFile, mode='w') as fp_out:
        result = hicstraw.straw(datatype, Norm, HiCFile, CHR1, CHR2, 'BP', resolution)
        for i in range(len(result)):
            print("{0}\t{1}\t{2}\t{3}\t{4}".format(chr1name, (result[i].binX + int(resolution / 2)), chr2name, (result[i].binY + int(resolution / 2)), result[i].counts), file=fp_out)
    
#===============================================
if __name__ == "__main__":
    main()

