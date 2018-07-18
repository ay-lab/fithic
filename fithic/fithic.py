#!/usr/bin/env python
'''
Created on Jan 29, 2014 by Ferhat Ay
Modified by Arya Kaul 2017-Present

'''

import sys
import math
import time
import numpy as np
from scipy import *
from scipy.interpolate import Rbf, UnivariateSpline
from scipy import optimize
import scipy.special as scsp
import bisect
import gzip
from scipy.stats.mstats import mquantiles
from scipy import stats
try:
    from . import myStats
    from . import myUtils
except:
    import myStats
    import myUtils
from sklearn.isotonic import IsotonicRegression
from sortedcontainers import SortedList
import os
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter, MaxNLocator
from pylab import *

versionStr = "You are using fithic version 2.0.3"

def parse_args(args):
    parser = argparse.ArgumentParser(description="Check the help flag")

    parser.add_argument("-i", "--interactions", dest="intersfile",\
                      help="REQUIRED: interactions between fragment pairs are \
                      read from INTERSFILE", required=True)

    parser.add_argument("-f", "--fragments", dest="fragsfile", \
                      help="REQUIRED: midpoints (or start indices) \
                      of the fragments are read from FRAGSFILE",\
                      required=True)

    parser.add_argument("-o", "--outdir", dest="outdir", \
                      help="REQUIRED: where the output files\
                      will be written", required=True)

    parser.add_argument("-r", "--resolution", dest="resolution", type=int,
                      help="REQUIRED: If the files are fixed size, please \
                      supply the resolution of the dataset here; otherwise, \
                      please use a value of 0 if the data is not fixed size." \
                      , required=True)
    
    parser.add_argument("-t", "--biases", dest="biasfile",\
                        help="RECOMMENDED: biases calculated by\
                        ICE or KR norm for each locus are read from BIASFILE",\
                        required=False)

    parser.add_argument("-p", "--passes", dest="noOfPasses",type=int,\
                        help="OPTIONAL: number of spline passes to run\
                        Default is 1",
                        required=False)

    parser.add_argument("-b", "--noOfBins", dest="noOfBins", type=int, \
                      help="OPTIONAL: number of equal-occupancy (count) \
                      bins. Default is 100", required=False)

    parser.add_argument("-m", "--mappabilityThres", dest="mappabilityThreshold",\
                      type=int, help="OPTIONAL: minimum number of hits per \
                      locus that has to exist to call it mappable. DEFAULT is 1.",\
                      required=False)

    parser.add_argument("-l", "--lib", dest="libname", help="OPTIONAL: Name of the\
                      library that is analyzed to be used for name of file prefixes \
                      . DEFAULT is fithic",
                      required=False)

    parser.add_argument("-U", "--upperbound", dest="distUpThres", type=int,
                      help="OPTIONAL: upper bound on the intra-chromosomal \
                      distance range (unit: base pairs). DEFAULT no limit. \
                      STRONGLY suggested to have a limit for large genomes,\
                      such as human/mouse. ex. '1000000, 5000000, etc.'",
                      required=False)

    parser.add_argument("-L", "--lowerbound", dest="distLowThres", type=int,
                      help="OPTIONAL: lower bound on the intra-chromosomal \
                      distance range (unit: base pairs). DEFAULT no limit. \
                      Suggested limit is 2x the resolution of the input files",
                      required=False)

    parser.add_argument("-v", "--visual", action="store_true", dest="visual",\
                      help="OPTIONAL: use this flag for generating plots. \
                      DEFAULT is False.", required=False)

    parser.add_argument("-x", "--contactType", dest="contactType",
                      help="OPTIONAL: use this flag to determine which chromosomal \
                      regions to study (intraOnly, interOnly, All) \
                      DEFAULT is intraOnly", required=False)
    
    parser.add_argument("-tL", "--biasLowerBound", dest="biasLowerBound",\
                      help="OPTIONAL: this flag is used to determine the lower bound\
                      of bias values to discard. DEFAULT is 0.5"\
                      , required=False)
    
    parser.add_argument("-tU", "--biasUpperBound", dest="biasUpperBound",\
                      help="OPTIONAL: this flag is used to determine the upper bound\
                      of bias values to discard. DEFAULT is 2"\
                      , required=False)

    parser.add_argument("-V", "--version", action="version",version=versionStr)

    return parser.parse_args()

def main():
    args = parse_args(sys.argv[1:])
    print("\n")
    print("GIVEN FIT-HI-C ARGUMENTS")
    print("=========================")

    ##PARSE REQUIRED ARGUMENTS##
    fragsFile = args.fragsfile
    if os.path.exists(fragsFile):
        print("Reading fragments file from: %s" % fragsFile)
    else:
        print("Fragment file not found")
        sys.exit(2)
    try:
        fragsF = gzip.open(fragsFile, 'r')
        fragsF.readline()
    except:
        print("Fragments file is not gzipped. Exiting now...")
        sys.exit(2)

    contactCountsFile = args.intersfile
    if os.path.isfile(contactCountsFile):
        print("Reading interactions file from: %s" % contactCountsFile)
    else:
        print("Interaction file not found")
        sys.exit(2)
    try:
        contactCountsF = gzip.open(contactCountsFile, 'r')
        contactCountsF.readline()
    except:
        print("Interactions file is not gzipped. Exiting now...")

    outputPath = args.outdir
    if not os.path.isdir(outputPath):
        os.makedirs(outputPath)
        print("Output path created %s" % outputPath)
    else:
        print("Output path being used from %s" % outputPath)

    resolution = args.resolution
    if args.resolution == 0:
        print("Fixed size data not being used.")
    elif args.resolution > 0:
        print("Fixed size option detected... Fast version of FitHiC will be used")
        print("Resolution is %s kb" % (resolution/1000))
    else:
        print("INVALID RESOLUTION ARGUMENT DETECTED")
        print("Please make sure the given resolution is a positive number greater than zero")
        print("User-given resolution: %s" % resolution)
        sys.exit(2)

    ##PARSE OPTIONAL ARGUMENTS##
    
    if args.biasfile is not None:
        if os.path.isfile(args.biasfile):
            print("Reading bias file from: %s" % args.biasfile)
        else:
            print("Bias file not found")
            sys.exit(2)
    else:
        print("No bias file")
    biasFile = args.biasfile 
    

    noOfPasses = 1
    if args.noOfPasses:
        noOfPasses = args.noOfPasses
    print("The number of spline passes is %s" % noOfPasses)

    noOfBins = 100
    if args.noOfBins:
        noOfBins = args.noOfBins
    print("The number of bins is %s" % noOfBins)

    global mappThres
    mappThres = 1
    if args.mappabilityThreshold:
        mappThres = args.mappabilityThreshold
    print("The number of reads required to consider an interaction is %s" % mappThres)

    libName = "FitHiC"
    if args.libname:
        libName = args.libname
    print("The name of the library for outputted files will be %s" % libName)

    global distLowThres
    global distUpThres
    distUpThres = float("inf")
    distLowThres = 0
    if args.distUpThres:
        distUpThres = args.distUpThres
    if args.distLowThres:
        distLowThres = args.distLowThres
    print("Upper Distance threshold is %s" % distUpThres)
    print("Lower Distance threshold is %s" % distLowThres)

    global visual
    visual = False
    if args.visual:
        #### matplotlib fontsize settings
        visual = True
        print("Graphs will be outputted")

    global interOnly
    global allReg
    chromosome_region=args.contactType
    if chromosome_region==None:
        chromosome_region='intraOnly'
    interOnly=False
    allReg=False
    if chromosome_region == "All":
        print("All genomic regions will be analyzed")
        allReg=True
    elif chromosome_region == "interOnly":
        print("Only inter-chromosomal regions will be analyzed")
        interOnly=True
    elif chromosome_region == "intraOnly":
        print("Only intra-chromosomal regions will be analyzed")
        interOnly=False
        allReg=False
    else:
        print("Invalid Option. Only options are 'All', 'interOnly', or 'intraOnly'")
        sys.exit(2)

    global biasLowerBound
    global biasUpperBound
    biasLowerBound = 0.5
    biasUpperBound = 2
    if args.biasLowerBound:
        biasLowerBound = args.biasLowerBound
    if args.biasUpperBound:
        biasUpperBound = args.biasUpperBound
    if biasLowerBound > biasUpperBound:
        print("Invalid Option. Bias lower bound is greater than bias upper bound. Please fix.")
        sys.exit(2)
    print("Lower bound of bias values is %s" % biasLowerBound)
    print("Upper bound of bias values is %s" % biasUpperBound)

    print("All arguments processed. Running FitHiC now...")
    print("=========================")
    print("\n")

    #########################PARSING COMPLETE############################################

    possibleIntraInRangeCount=0 # count of all possible in range intra-chr fragpairs
    observedIntraInRangeCount=0 # count of obs. in range intra-chr frags based on intxn file
    observedIntraInRangeSum=0 # sum of all observed intra-chr read counts in range

    possibleIntraAllCount=0 # Same as above, but without range restriction
    observedIntraAllCount=0
    observedIntraAllSum=0

    possibleInterAllCount=0 # Same as above, note that the notion of distance thresholds
                            # does not apply for interchr intxns
    observedInterAllCount=0
    observedInterAllSum=0

    global baselineIntraChrProb
    baselineIntraChrProb=0  # 1.0/possibleIntraAllCount
    global interChrProb 
    interChrProb=0 # 1.0/possibleInterAllCount

    minObservedGenomicDist=float("inf")
    maxObservedGenomicDist=0
    maxPossibleGenomicDist=0

    #distScaling just avoids overflow - but is necessary for large genomes
    global distScaling
    distScaling=1000000.0 #RUNBY

    global toKb
    global toMb
    global toProb
    toKb=10**-3
    toMb=10**-6
    toProb=10**5

    #intermediate values outputted here
    global logfile
    logfile = os.path.join(outputPath, libName+".fithic.log")

    ##maindic will be generated first using the interactions file only
    mainDic={} # given a distance this dictionary will return [Npairs,TotalContactCount] for only those interactions present in the interactions file
    (mainDic,observedInterAllSum,observedIntraAllSum,observedIntraInRangeSum) = read_Interactions(contactCountsFile, biasFile)
    binStats = makeBinsFromInteractions(mainDic, noOfBins, observedIntraInRangeSum)

    #Enumerate (fast version) or generate (otherwise) all possible pairs of fragments within the range of interest.
    (binStats,noOfFrags, maxPossibleGenomicDist, possibleIntraInRangeCount, possibleInterAllCount, interChrProb, baselineIntraChrProb)= generate_FragPairs(binStats, fragsFile, resolution)

    #read and parse bias values for each locus from ICE or KR normalization output
    if biasFile:
        biasDic = read_biases(biasFile)
    else:
        biasDic = 0

    #bin the data in desired number of bins, and for each bin, calculate the average genomic distance and average contact probability
    (x,y,yerr)= calculateProbabilities(mainDic, binStats,resolution,os.path.join(outputPath,libName+".fithic_pass1"), observedIntraInRangeSum)

    splinefit1st=time.time()
    print("Spline fit Pass 1 starting...")
    outliersline = SortedList()
    outliersdist = SortedList()
    #fit a smooth spline to the bin values, and compute and write p values/q values
    splineXinit,splineYinit,residual,outliersline, outliersdist, FDRXinit, FDRYinit= fit_Spline(mainDic,x,y,yerr,contactCountsFile,os.path.join(outputPath,libName+".spline_pass1"),biasDic, outliersline, outliersdist, observedIntraInRangeSum, possibleIntraInRangeCount, possibleInterAllCount, observedIntraAllSum, observedInterAllSum, resolution, 1)
    print("Number of outliers is... %s" % len(outliersline))
    splinefit1en = time.time()
    print("Spline fit Pass 1 completed. Time took %s" % (splinefit1en-splinefit1st))

    ### DO THE NEXT PASSES IF REQUESTED ###
    for i in range(2,1+noOfPasses):
        if interOnly:
            print("Extra spline fits will not help with interOnly spline fit... Bypassing option")
            break
        print("\n")
        print("\n")
        (mainDic,observedInterAllSum,observedIntraAllSum,observedIntraInRangeSum) = read_Interactions(contactCountsFile, biasFile, outliersline)
        binStats = makeBinsFromInteractions(mainDic, noOfBins, observedIntraInRangeSum, outliersdist)
        (binStats,noOfFrags, maxPossibleGenomicDist, possibleIntraInRangeCount,possibleInterAllCount, interChrProb, baselineIntraChrProb)= generate_FragPairs(binStats, fragsFile, resolution)
        (x,y,yerr)= calculateProbabilities(mainDic, binStats,resolution,os.path.join(outputPath,libName+".fithic_pass"+str(i)), observedIntraInRangeSum)
        splinefitst=time.time()
        print("Spline fit Pass %s starting..." % i)
        splineX,splineY,residual,outliersline, outliersdist, FDRX, FDRY= fit_Spline(mainDic,x,y,yerr,contactCountsFile,os.path.join(outputPath,libName+".spline_pass"+str(i)),biasDic, outliersline, outliersdist, observedIntraInRangeSum, possibleIntraInRangeCount, possibleInterAllCount, observedIntraAllSum, observedInterAllSum, resolution, i)
        splinefiten = time.time()
        print("Spline fit Pass %s completed. Time took %s" % (i,(splinefit1en-splinefit1st)))
        if visual:
            compare_Spline_FDR(FDRXinit, FDRYinit, FDRX, FDRY, os.path.join(outputPath, libName+".spline_FDR_comparison"),str(i))
            compareFits_Spline(splineXinit, splineYinit, splineX, splineY, os.path.join(outputPath,libName+".spline_comparison"), str(i))
    print("=========================")
    print("Fit-Hi-C completed successfully")
    print("\n")

##FUNCTIONS START###

def read_Interactions(contactCountsFile, biasFile, outliers=None):
    mainDic={}
    print("Reading the contact counts file to generate bins...")
    startT = time.time()

    observedInterAllSum=0 #used
    observedIntraAllSum=0 #used
    observedInterAllCount=0
    observedIntraAllCount=0 #notused
    observedIntraInRangeSum=0 #used
    observedIntraInRangeCount=0 #notused
    minObservedGenomicDist=float('inf') #notused
    maxObservedGenomicDist=0 #notused

    linectr = 0
    outlierposctr = 0
    #Loop through every line in the contactCountsFile
    with gzip.open(contactCountsFile, 'rt') as f:
        for lines in f:
            if outliers != None and outlierposctr<len(outliers):
                if linectr == outliers[outlierposctr]:
                    linectr+=1
                    outlierposctr+=1
                    continue
            ch1,mid1,ch2,mid2,contactCount=lines.split()
            #create the interaction
            contactCount=float(contactCount)
            interxn=myUtils.Interaction([ch1, int(mid1), ch2, int(mid2)])
            interxn.setCount(contactCount)
            interactionType = interxn.getType(distLowThres,distUpThres)
            if interactionType=='inter':
                observedInterAllSum += interxn.getCount()
                observedInterAllCount +=1
            else: # any type of intra
                observedIntraAllSum +=interxn.getCount()
                observedIntraAllCount +=1
                if interactionType=='intraInRange':
                    #interxn.setDistance(interxn.getDistance()+(1000-interxn.getDistance()) % 1000)
                    minObservedGenomicDist=min(minObservedGenomicDist,interxn.getDistance())
                    maxObservedGenomicDist=max(maxObservedGenomicDist,interxn.getDistance())
                    if interxn.getDistance() not in mainDic:
                        mainDic[interxn.getDistance()] = [0,0]
                    mainDic[interxn.getDistance()][1]+=interxn.getCount()
                    observedIntraInRangeSum +=interxn.getCount()
                    observedIntraInRangeCount +=1
            linectr+=1
    endT = time.time()
    print("Interactions file read. Time took %s" % (endT-startT))
    with open(logfile, 'w') as log:
        log.write("\n\nInteractions file read successfully\n")
        log.write("------------------------------------------------------------------------------------\n")
        log.write("Observed, Intra-chr in range: pairs= "+str(observedIntraInRangeCount) +"\t totalCount= "+str(observedIntraInRangeSum)+"\n")
        log.write("Observed, Intra-chr all: pairs= "+str(observedIntraAllCount) +"\t totalCount= "+str(observedIntraAllSum)+"\n")
        log.write("Observed, Inter-chr all: pairs= "+str(observedInterAllCount) +"\t totalCount= "+str(observedInterAllSum)+"\n")
        log.write("Range of observed genomic distances [%s %s]" % (minObservedGenomicDist,maxObservedGenomicDist) + "\n"),
        log.write("\n")
    return (mainDic,observedInterAllSum,observedIntraAllSum,observedIntraInRangeSum) # from read_Interactions

def makeBinsFromInteractions(mainDic,noOfBins, observedIntraInRangeSum, outliersdist=None):
    with open(logfile, 'a') as log:
        log.write("Making equal occupancy bins\n")
        log.write("------------------------------------------------------------------------------------\n")
        noPerBin = observedIntraInRangeSum/noOfBins
        log.write("Observed intra-chr read counts in range\t"+repr(observedIntraInRangeSum)+ "\nDesired number of contacts per bin\t" +repr(noPerBin)+",\nNumber of bins\t"+repr(noOfBins)+"\n")

    # the following five lists will be the print outputs
    interactionTotalForBinTermination=0
    n=0 # bin counter so far
    totalInteractionCountSoFar=0
    distsToGoInAbin=[]
    binFull=0
    desiredPerBin=(observedIntraInRangeSum)/noOfBins
    bins = []
    for i in sorted(mainDic.keys()): #everything here is inrange by definition
        totalInteractionCountSoFar+=mainDic[i][1]
        # if one distance has more than necessary counts to fill a bin
        if mainDic[i][1]>=desiredPerBin:
            distsToGoInAbin.append(i)
            interactionTotalForBinTermination=0
            binFull=1
        # if adding the next bin will fill the bin
        elif interactionTotalForBinTermination+mainDic[i][1] >= desiredPerBin:
            distsToGoInAbin.append(i)
            interactionTotalForBinTermination=0
            binFull=1
        # if adding the next bin will not fill the bin
        else:
            distsToGoInAbin.append(i)
            interactionTotalForBinTermination+=mainDic[i][1]
        # if bin is already full
        if binFull==1:
            noOfPairsForBin=0
            interactionTotalForBin=0
            avgDistance=0
            # dynamically update the desiredPerBin after each bin is full
            n+=1
            if n<noOfBins:
                desiredPerBin=1.0*(observedIntraInRangeSum-totalInteractionCountSoFar)/(noOfBins-n)
            bins.append(distsToGoInAbin)
            interactionTotalForBinTermination=0
            binFull=0
            distsToGoInAbin=[]
    #print(bins)
    binStats = {}
    for binIdx in range(len(bins)):
        ##binStats
        #0: range of distances in this bin
        #1: no. of possible pairs w/in this range of distances
        #2: sumoverallContactCounts
        #3: Sumoveralldistances in this bin in distScaling vals
        #4: avg CC
        #5: avg distance
        #6: bins
        if binIdx == 0:
            lb = 0
        else:
            lb = max(bins[binIdx-1])+1
        ub = bins[binIdx][-1]
        binStats[binIdx]=[(lb, ub), 0, 0, 0, 0, 0, bins[binIdx], 0]
        for dists in bins[binIdx]:
            binStats[binIdx][2]+=mainDic[dists][1]
            #binStats[binIdx][3]+=(dists/distScaling)
    
    if outliersdist != None:
        binTracker = 0
        for i in range(len(outliersdist)):
            intxnDistance = outliersdist[i]
            currBin = binStats[binTracker]
            minOfBin = currBin[0][0]
            maxOfBin = currBin[0][1]
            while not (minOfBin<=intxnDistance<=maxOfBin):
                binTracker += 1
                if binTracker not in binStats:
                    binTracker-=1
                    currBin = binStats[binTracker]
                    minOfBin = currBin[0][0]
                    maxOfBin = currBin[0][1]
                    break
                else:
                    currBin = binStats[binTracker]
                    minOfBin = currBin[0][0]
                    maxOfBin = currBin[0][1]
            currBin[7]-=1
            currBin[1]-=1
    
    with open(logfile, 'a') as log:
        log.write("Equal occupancy bins generated\n")
        log.write("\n")
    return binStats

def generate_FragPairs(binStats, fragsfile, resolution): 
    if resolution:
        with open(logfile, 'a') as log:
            log.write("Looping through all possible fragment pairs in-range\n")
            log.write("------------------------------------------------------------------------------------\n"),
    else:
        with open(logfile, 'a') as log:
            log.write("Enumerating all possible fragment pairs in-range\n")
            log.write("------------------------------------------------------------------------------------\n"),
    startT = time.time()

    minPossibleGenomicDist = float("inf")
    maxPossibleGenomicDist = 0
    possibleIntraAllCount = 0
    possibleInterAllCount = 0
    possibleIntraInRangeCount = 0
    interChrProb = 0
    baselineIntraChrProb = 0

    allFragsDic={}
    with gzip.open(fragsfile,'rt') as infile:
        for line in infile:
            words=line.split()
            currChr=words[0]
            currMid=int(words[2])
            currHit=int(words[3])
            if currChr not in allFragsDic:
                allFragsDic[currChr]=[]
            if currHit>=mappThres:
                allFragsDic[currChr].append(currMid)

    if resolution:
        noOfFrags=0
        maxFrags={}
        
        for ch in allFragsDic:
            maxFrags[ch]=max([int(i)-resolution/2 for i in allFragsDic[ch]])
            noOfFrags+=len(allFragsDic[ch])
            maxPossibleGenomicDist=max(maxPossibleGenomicDist,maxFrags[ch])
        
        for ch in sorted(allFragsDic.keys()):
            maxFrag=maxFrags[ch]
            n=len(allFragsDic[ch])
            d=0
            binTracker = 0
            possibleIntraInRangeCountPerChr = 0
            for intxnDistance in range(0,int(maxFrag+1),resolution):
                npairs = n-d
                d+=1
                if myUtils.in_range_check(intxnDistance,distLowThres,distUpThres):
                    minPossibleGenomicDist = min(minPossibleGenomicDist, intxnDistance)
                    possibleIntraInRangeCountPerChr += npairs
                else:
                    continue
                currBin = binStats[binTracker]
                minOfBin = currBin[0][0]
                maxOfBin = currBin[0][1]
                while not (minOfBin<=intxnDistance<=maxOfBin):
                    binTracker += 1
                    if binTracker not in binStats:
                        binTracker-=1
                        currBin = binStats[binTracker]
                        minOfBin = currBin[0][0]
                        maxOfBin = currBin[0][1]
                        break
                    else:
                        currBin = binStats[binTracker]
                        minOfBin = currBin[0][0]
                        maxOfBin = currBin[0][1]
                currBin[7]+=npairs
                currBin[1]+=npairs
                currBin[3]+=(float(intxnDistance/distScaling)*npairs)
               # possibleIntraInRangeCountPerChr += npairs
            possibleInterAllCount+=n*(noOfFrags-n)
            possibleIntraAllCount+=(n*(n+1))/2 # n(n-1) if excluding self
            with open(logfile, 'a') as log:
                log.write("Chromosome " +repr(ch) +",\t"+str(n) +" mappable fragments, \t"+str(possibleIntraInRangeCountPerChr)\
                +" possible intra-chr fragment pairs in range,\t" + str((noOfFrags-n)*n) +" possible inter-chr fragment pairs\n")
            possibleIntraInRangeCount += possibleIntraInRangeCountPerChr
        possibleInterAllCount/=2
        try:
            interChrProb=1.0/possibleInterAllCount
        except:
            interChrProb = 0
        baselineIntraChrProb=1.0/possibleIntraAllCount

    else:
        noOfFrags = 0
        for ch in allFragsDic:
            noOfFrags += len(allFragsDic[ch])

        for ch in sorted(allFragsDic.keys()):
            countIntraPairs = 0
            fragsPerChr = sorted(allFragsDic[ch])
            templen = len(fragsPerChr)
            possibleInterAllCount += (noOfFrags-templen)*templen
            possibleIntraInRangeCountPerChr = 0
            for x in range(templen):
                binTracker = 0
                d = 0
                for y in range(x+1,templen):
                    intxnDistance = abs(float(fragsPerChr[x])-float(fragsPerChr[y]))
                    if myUtils.in_range_check(intxnDistance, distLowThres,distUpThres):
                        possibleIntraInRangeCountPerChr += 1 
                    else:
                        continue
                    maxPossibleGenomicDist = max(maxPossibleGenomicDist, intxnDistance)
                    minPossibleGenomicDist = min(minPossibleGenomicDist, intxnDistance)
                    npairs = templen-d
                    d+=1
                    currBin = binStats[binTracker]
                    minOfBin = currBin[0][0]
                    maxOfBin = currBin[0][1]
                    while not (minOfBin<=intxnDistance<=maxOfBin):
                        binTracker += 1
                        if binTracker not in binStats:
                            binTracker-=1
                            currBin = binStats[binTracker]
                            minOfBin = currBin[0][0]
                            maxOfBin = currBin[0][1]
                            break
                        else:
                            currBin = binStats[binTracker]
                            minOfBin = currBin[0][0]
                            maxOfBin = currBin[0][1]
                    currBin[7]+=npairs
                    currBin[1]+=1
                    currBin[3]+=float(intxnDistance/distScaling)*npairs
                    possibleIntraAllCount += 1
            with open(logfile, 'a') as log:
                log.write("Chromosome " +repr(ch) +",\t"+str(templen) +" mappable fragments, \t"+str(possibleIntraInRangeCountPerChr)\
                +" possible intra-chr fragment pairs in range,\t" + str((noOfFrags-templen)*templen) +" possible inter-chr fragment pairs\n")
            possibleIntraInRangeCount += possibleIntraInRangeCountPerChr
        possibleInterAllCount/=2
        try:
            interChrProb=1.0/possibleInterAllCount
        except:
            interChrProb = 0
        baselineIntraChrProb=1.0/possibleIntraAllCount
    endT = time.time()
    print("Fragments file read. Time took %s" % (endT-startT))

    with open(logfile, 'a') as log:
        log.write("Number of all fragments= %s\n" % (noOfFrags))
        log.write("Possible, Intra-chr in range: pairs= %s \n" % (possibleIntraInRangeCount))
        log.write("Possible, Intra-chr all: pairs= %s \n" % (possibleIntraAllCount))
        log.write("Possible, Inter-chr all: pairs= %s \n" % (possibleInterAllCount))
        log.write("Desired genomic distance range   [%d %s] \n" % (distLowThres,distUpThres)),
        log.write("Range of possible genomic distances  [%d  %d] \n" % (minPossibleGenomicDist, maxPossibleGenomicDist)),
        log.write("Baseline intrachromosomal probability is %s \n" % (baselineIntraChrProb)),
        log.write("Interchromosomal probability is %s \n" % (interChrProb)),

    return (binStats,noOfFrags, maxPossibleGenomicDist, possibleIntraInRangeCount, possibleInterAllCount, interChrProb, baselineIntraChrProb) # return from generate_FragPairs


def read_biases(infilename):
    startt = time.time()
    biasDic={}

    rawBiases=[]
    try:
        infile =gzip.open(infilename, 'rt')
    except:
        infile = open(infilename, 'r')
    for line in infile:
        words=line.rstrip().split()
        chr=words[0]; midPoint=int(words[1]); bias=float(words[2])
        if bias!=1.0:
           rawBiases.append(bias)
    botQ,med,topQ=mquantiles(rawBiases,prob=[0.05,0.5,0.95])
    with open(logfile, 'a') as log:
        log.write("5th quantile of biases: "+str(botQ)+"\n")
        log.write("50th quantile of biases: "+str(med)+"\n")
        log.write("95th quantile of biases: "+str(topQ)+"\n")
    totalC=0
    discardC=0
    for line in infile:
        words=line.rstrip().split()
        chr=words[0]; midPoint=int(words[1]); bias=float(words[2]);
        if bias<biasLowerBound:
            bias=-1 #botQ
            discardC+=1
        elif bias>biasUpperBound:
            bias=-1 #topQ
            #bias=1
            discardC+=1
        totalC+=1
        if chr not in biasDic:
            biasDic[chr]={}
        if midPoint not in biasDic[chr]:
            biasDic[chr][midPoint]=bias
    infile.close()
    with open(logfile, 'a') as log:
        log.write("Out of " + str(totalC) + " loci " +str(discardC) +" were discarded with biases not in range [0.5 2]\n\n" )
    endt = time.time()
    print("Bias file read. Time took %s" % (endt-startt))
    return biasDic # from read_biases

def calculateProbabilities(mainDic,binStats,resolution,outfilename,observedIntraInRangeSum):
    with open(logfile, 'a') as log:
        log.write("\nCalculating probability means and standard deviations of contact counts\n"),
        log.write("------------------------------------------------------------------------------------\n"),

    if resolution:
        nameoffile = (outfilename+'.res'+str(resolution)+'.txt')
    else:
        nameoffile = (outfilename+'.txt')

    outfile=open(nameoffile, 'w')
    x = []
    y = []
    yerr = []
    pairCounts=[]
    interactionTotals=[]
        ##binStats
        #0: range of distances in this bin
        #1: no. of possible pairs w/in this range of distances
        #2: sumoverallContactCounts
        #3: Sumoveralldistances in this bin in distScaling vals
        #4: avg CC
        #5: avg distance
        #6: bins
        #7: no. of possible pairs w/ proper dist

    for i in range(len(binStats)):
        currBin = binStats[i]
        sumCC = currBin[2]
        sumDistB4Scaling = currBin[3]
        possPairsInRange = currBin[1]
        try:
            avgCC = (1.0*sumCC/possPairsInRange)/observedIntraInRangeSum
        except:
            print("WARNING - Zero avg. contact in bin. Ensure interaction file is correct.")
            avgCC = 0
        try:
            avgDist = distScaling*(sumDistB4Scaling/currBin[7])
        except:
            print("WARNING - Zero avg. distance in bin. Ensure interaction file is correct.")
            avgDist = 0
        currBin[4]=avgCC
        currBin[5]=avgDist
        y.append(avgCC)
        x.append(avgDist)
        """
        meanCountPerPair = 0
        M2 = 0
        for dists in currBin[6]: #by definition not including the nonzero dists in this bin in this calc.
            delta = mainDic[dists][1]-meanCountPerPair
            meanCountPerPair += (delta*1.0)/possPairsInRange
            M2 += delta*(mainDic[dists][1]-meanCountPerPair)
        var = M2/(possPairsInRange-1)
        sd = math.sqrt(var)
        se = sd/math.sqrt(possPairsInRange)
        se_p = se/observedIntraInRangeSum
        #yerr.append(se_p)
        """
        yerr.append(0)
        pairCounts.append(possPairsInRange)
        interactionTotals.append(sumCC)

    print("Writing %s" % nameoffile)
    outfile.write("avgGenomicDist\tcontactProbability\tstandardError\tnoOfLocusPairs\ttotalOfContactCounts\n")
    for i in range(len(x)):
        outfile.write("%d" % x[i] + "\t"+"%.2e" % y[i]+ "\t" + "%.2e" % yerr[i] + "\t" +"%d" % pairCounts[i] + "\t" +"%d" % interactionTotals[i]+"\n")
    outfile.close()
    with open(logfile, 'a') as log:
        log.write("Means and error written to %s\n" % (nameoffile)),
        log.write("\n"),
    return [x,y,yerr] # from calculateProbabilities


def fit_Spline(mainDic,x,y,yerr,infilename,outfilename,biasDic,outliersline,outliersdist,observedIntraInRangeSum, possibleIntraInRangeCount, possibleInterAllCount, observedIntraAllSum, observedInterAllSum, resolution, passNo):
    with open(logfile, 'a') as log:
        log.write("\nFitting a univariate spline to the probability means\n"),
        log.write("------------------------------------------------------------------------------------\n"),
   
    splineX = None
    newSplineY = None
    residual = None 
    FDRx = None
    FDRy = None

    if not interOnly:
        if outliersdist != None:
            y = [f for _, f in sorted(zip(x,y), key=lambda pair: pair[0])]
            x.sort()
        for i in range(1,len(x)):
            if x[i]<=x[i-1]:
                print("ERROR in spline fitting. Distances do not decrease across bins. Ensure interaction file is correct.")
                print("Avg. distance of bin(i-1)... %s" % x[i-1])
                print("Avg. distance of bin(i)... %s" % x[i])
                sys.exit(2)
        
        # maximum residual allowed for spline is set to min(y)^2
        splineError=min(y)*min(y)

        # use fitpack2 method -fit on the real x and y from equal occupancy binning
        ius = UnivariateSpline(x, y, s=splineError)
        tempMaxX=max(x)
        tempMinX=min(x)
        tempList=sorted([dis for dis in mainDic])
        splineX=[]
        ### The below for loop will make sure nothing is out of range of [min(x) max(x)]
        ### Therefore everything will be within the range where the spline is defined
        for i in tempList:
            if tempMinX<=i<=tempMaxX:
                splineX.append(i)
        splineY=ius(splineX)
        #print(splineY)
        #print(yerr)


        ir = IsotonicRegression(increasing=False)
        newSplineY = ir.fit_transform(splineX,splineY)
        #print(newSplineY)
        residual =sum([i*i for i in (y - ius(x))])

        if visual==True:
            xi = np.linspace(min(x),max(x),5*len(x))
            yi = ius(xi)

            print("Plotting %s" % (outfilename + ".png"))
            plt.clf()
            fig = plt.figure()
            ax = fig.add_subplot(2,1,1)
            plt.plot(myUtils.scale_a_list(splineX,toKb), myUtils.scale_a_list(newSplineY,toProb),'g-',label="spline-"+str(passNo),linewidth=2)
            plt.errorbar(myUtils.scale_a_list(x,toKb),myUtils.scale_a_list(y,toProb),myUtils.scale_a_list(yerr,toProb),fmt='r.',label="Mean with std. error",linewidth=2) 
        
            #plt.ylabel('Contact probability (x10$^{-5}$)',fontsize='large')
            #plt.xlabel('Genomic distance (kb)',fontsize='large')
            plt.ylabel('Contact probability (x10$^{-5}$)')
            plt.xlabel('Genomic distance (kb)')
            if distLowThres>0 and distUpThres<float("inf"):
                plt.xlim(myUtils.scale_a_list([distLowThres, distUpThres],toKb))
            plt.gca().yaxis.set_major_locator( MaxNLocator(nbins = 3, prune=None))
            ax.legend(loc="upper right")

            ax = fig.add_subplot(2,1,2)

            plt.loglog(splineX,newSplineY,'g-')
            plt.errorbar(x, y, yerr=yerr, fmt='r.') # Data
            if distLowThres>0 and distUpThres<float("inf"):
                plt.xlim([distLowThres, distUpThres])
            plt.ylabel('Contact probability (log-scale)')
            plt.xlabel('Genomic distance (log-scale)')

            plt.savefig(outfilename+'.png')
            

    # NOW write the calculated pvalues and corrected pvalues in a file
    infile = gzip.open(infilename, 'rt')
    intraInRangeCount=0
    intraOutOfRangeCount=0
    intraVeryProximalCount=0
    interCount=0
    discardCount=0
    p_vals=[]
    q_vals=[]
    biasl=[]
    biasr=[]
    for line in infile:
        ch1,mid1,ch2,mid2,contactCount=line.rstrip().split()
        contactCount = float(contactCount)
        interxn=myUtils.Interaction([ch1, int(mid1), ch2, int(mid2)])
        interxn.setCount(contactCount)
        mid1 = int(mid1); mid2 = int(mid2)
        interactionType = interxn.getType(distLowThres,distUpThres)
        bias1=1.0; bias2=1.0;  # assumes there is no bias to begin with
        # if the biasDic is not null sets the real bias values
        if biasDic:
            if ch1 in biasDic and mid1 in biasDic[ch1]:
                bias1=biasDic[ch1][mid1]
            if ch2 in biasDic and mid2 in biasDic[ch2]:
                bias2=biasDic[ch2][mid2]
        biasl.append(bias1)
        biasr.append(bias2)
        if (bias1<0 or bias2<0) and interactionType !='inter':
            prior_p=1.0
            p_val=1.0
            discardCount+=1
        elif interactionType=='intraInRange' and not interOnly:
            distToLookUp=max(interxn.getDistance(),min(x))
            distToLookUp=min(distToLookUp,max(x))
            i=min(bisect.bisect_left(splineX, distToLookUp),len(splineX)-1)
            prior_p=newSplineY[i]*(bias1*bias2) 
            p_val=scsp.bdtrc(interxn.getCount()-1,observedIntraInRangeSum,prior_p)
            intraInRangeCount +=1
        elif interactionType =='intraShort' and not interOnly:
            prior_p=1.0
            p_val=1.0
            intraVeryProximalCount += 1
        elif interactionType =='intraLong' and not interOnly:
            prior_p=1.0
            #p_val=scsp.bdtrc(interxn.getCount()-1, observedIntraAllSum,prior_p) ##RUNBY
            p_val=1.0
            intraOutOfRangeCount += 1
        else:
            if allReg or interOnly:
                prior_p=interChrProb*(bias1*bias2)
                p_val=scsp.bdtrc(interxn.getCount()-1,observedInterAllSum,prior_p)
                interCount += 1
            else:
                p_val=1.0
                #p_vals.append(p_val)
        p_vals.append(p_val)
    infile.close()

    outlierThres = 0
    # Do the BH FDR correction
    if allReg:
        outlierThres=1.0/(possibleIntraInRangeCount+possibleInterAllCount)
        q_vals=myStats.benjamini_hochberg_correction(p_vals, possibleInterAllCount+possibleIntraInRangeCount)
    elif interOnly and not allReg:
        outlierThres = 1.0/possibleInterAllCount
        q_vals=myStats.benjamini_hochberg_correction(p_vals, possibleInterAllCount)
    else:
        outlierThres = 1.0/possibleIntraInRangeCount
        q_vals=myStats.benjamini_hochberg_correction(p_vals, possibleIntraInRangeCount)
    print("Outlier threshold is... %s" % (outlierThres))

    #now we write the values back to the file
    infile =gzip.open(infilename, 'rt')
    if resolution:
        outfile =gzip.open(outfilename+'.res'+str(resolution)+'.significances.txt.gz', 'wt')
    else:
        outfile =gzip.open(outfilename+'.significances.txt.gz', 'wt')
    print("Writing p-values and q-values to file %s" % (outfilename + ".significances.txt"))
    outfile.write("chr1\tfragmentMid1\tchr2\tfragmentMid2\tcontactCount\tp-value\tq-value\tbias1\tbias2\n")
    count=0
    for line in infile:
        words=line.rstrip().split()
        chr1=words[0]
        midPoint1=int(words[1])
        chr2=words[2]
        midPoint2=int(words[3])
        interactionCount=float(words[4])
        p_val=p_vals[count]
        q_val=q_vals[count]
        bias1=biasl[count]
        bias2=biasr[count]
        
        if (allReg or interOnly) and chr1!=chr2:
            outfile.write("%s\t%d\t%s\t%d\t%d\t%e\t%e\t%e\t%e\n" % (str(chr1), midPoint1, str(chr2), midPoint2, interactionCount, p_val, q_val, bias1, bias2))
        if (allReg or not interOnly) and chr1==chr2:
            interactionDistance = abs(midPoint1-midPoint2)
            if myUtils.in_range_check(interactionDistance,distLowThres, distUpThres):
                outfile.write("%s\t%d\t%s\t%d\t%d\t%e\t%e\t%e\t%e\n" % (str(chr1), midPoint1, str(chr2), midPoint2, interactionCount, p_val, q_val, bias1, bias2))
        
        if p_val<outlierThres:
            outliersline.add(count)
            outliersdist.add(abs(midPoint1-midPoint2))
        count+=1
    outfile.close()
    infile.close()
    if visual == True:
        print("Plotting q-values to file %s" % outfilename + ".qplot.png")
    minFDR=0.0
    maxFDR=0.05
    increment=0.001
    FDRx,FDRy=plot_qvalues(q_vals,minFDR,maxFDR,increment,outfilename+".qplot")
        
    with open(logfile, 'a') as log:
        log.write("Spline successfully fit\n"),
        log.write("\n"),
        log.write("\n"),

    return [splineX, newSplineY, residual, outliersline, outliersdist, FDRx, FDRy] # from fit_Spline

def plot_qvalues(q_values,minFDR,maxFDR,increment,outfilename):
    qvalTicks=np.arange(minFDR,maxFDR+increment,increment)
    significantTicks=[0 for i in range(len(qvalTicks))]
    qvalBins=[-1 for i in range(len(q_values))]
    for i, q in enumerate(q_values):
        if math.isnan(q): q=1 #make sure NaNs are set to 1 
        qvalBins[i]=int(math.floor(q/increment))
    
    for i in range(len(qvalBins)):
        if qvalBins[i]>=len(qvalTicks):
            continue
        significantTicks[qvalBins[i]]+=1
    
    # make it cumulative 
    for i in range(1,len(significantTicks)):
        significantTicks[i]=significantTicks[i]+significantTicks[i-1]
    # shift them by 1
    for i in range(1,len(significantTicks)):
        significantTicks[-1*i]=significantTicks[-1*i-1]
    significantTicks[0]=0

    if visual==True:
        plt.clf()
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        plt.plot(qvalTicks,significantTicks, 'b*-')
        plt.xlabel('FDR threshold')
        plt.ylabel('Number of significant contacts')
        plt.savefig(outfilename+'.png')

    return [qvalTicks,significantTicks]

def compare_Spline_FDR(splineFDRxinit,splineFDRyinit,splineFDRx,splineFDRy,figname,i):
    newlab = 'spline-' + str(i)
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    plt.plot(splineFDRx[1:],myUtils.scale_a_list(splineFDRy[1:],toKb), 'r+-',label=newlab)
    plt.plot(splineFDRxinit[1:],myUtils.scale_a_list(splineFDRyinit[1:],toKb), 'g.-',label='spline-1')
    plt.xlabel('FDR threshold')
    plt.ylabel('Significant contacts (x10$^{3}$)')
    plt.gca().yaxis.set_major_locator( MaxNLocator(prune='lower'))
    lg=ax.legend(loc="lower right")
    lg.draw_frame(False)
    plt.savefig(figname+'.png')

def compareFits_Spline(splineXinit,splineYinit,splineX,splineY,figname,X):
    downsample=min(5000,len(splineXinit))
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    x=splineXinit
    y=splineYinit
    col='g.-'
    lab='spline-1'
    randIndcs=[i for i in range(len(x))]
    randIndcs=np.random.choice(randIndcs,downsample)
    randIndcs=sorted(randIndcs)
    x=myUtils.scale_a_list([x[i] for i in randIndcs],toKb)
    y=myUtils.scale_a_list([y[i] for i in randIndcs],toProb)
    plt.plot(x,y,col,label=lab)

    if figname[-1]!='1': # meaning this is not the very first step
        x=splineX
        y=splineY
        col='r.-'
        lab='spline-'+X
        randIndcs=[i for i in range(len(x))]
        randIndcs=np.random.choice(randIndcs,downsample)
        randIndcs=sorted(randIndcs)
        x=myUtils.scale_a_list([x[i] for i in randIndcs],toKb)
        y=myUtils.scale_a_list([y[i] for i in randIndcs],toProb)
        plt.plot(x,y,col,label=lab)
    else: # plot only at a limited range and plot discrete binning
        if max(x)>1000: # if it's a big genome
            plt.xlim([500,1000])
            plt.ylim([0,1.0])
        else: # small genome
            plt.xlim([50,100])
            plt.ylim([0,0.5])

    ax.legend(loc="upper right")
    plt.xlabel('Genomic distance (kb)')
    plt.ylabel('Contact probability (x10$^{-5}$)')
    plt.gca().yaxis.set_major_locator( MaxNLocator(prune='lower'))
    plt.savefig(figname+'.png')

if __name__ == "__main__":
    main()
