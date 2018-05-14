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
import myStats
import myUtils
from sklearn.isotonic import IsotonicRegression
import os
import argparse

versionStr = "You are using fithic version 2.0.0"

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
    
    parser.add_argument("-t", "--biases", dest="biasfile",\
                        help="OPTIONAL: biases calculated by\
                        ICE or KR norm for each locus are read from BIASFILE",\
                        required=False)
    
    parser.add_argument("-p", "--passes", dest="noOfPasses",type=int,\
                        help="OPTIONAL: number of passes after the initial\
                        (before) fit. If 'fast' default is 0, otherwise 1", 
                        required=False)

    parser.add_argument("-b", "--noOfBins", dest="noOfBins", type=int, \
                      help="OPTIONAL: number of equal-occupancy (count) \
                      bins. Default is 100", required=False)

    parser.add_argument("-m", "--mappabilityThres", dest="mappabilityThreshold",\
                      type=int, help="OPTIONAL: minimum number of hits per \
                      locus that has to exist to call it mappable. DEFAULT is 1.",\
                      required=False)
    
    parser.add_argument("-l", "--lib", dest="libname", help="OPTIONAL: Name of the\
                      library that is analyzed to be used for plots. DEFAULT is fithic",
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
    
    parser.add_argument("-x", "--chromosome_region", dest="chromosome_region", 
                      help="OPTIONAL: use this flag to determine which chromosomal \
                      regions to study (intraOnly, interOnly, All) \
                      DEFAULT is intraOnly", required=False)
    
    parser.add_argument("-r", "--resolution", dest="resolution", type=int,  
                      help="OPTIONAL: Use this option if the fragments file \
                      is fixed size binned data). DEFAULT is None", required=False)
    
    parser.add_argument("-V", "--version", action="version",version=versionStr)
    
    return parser.parse_args()

def main():
    args = parse_args(sys.argv[1:])

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
   
    ##PARSE OPTIONAL ARGUMENTS##
    biasFile = None
    if args.biasfile:
        if os.path.isfile(args.biasfile):
            print("Reading bias file from: %s" % args.biasfile)
        else:
            print("Bias file not found")
            sys.exit(2)
        biasFile = args.biasfile
    else:
        print("No bias file being used")
    
    noOfPasses = 0
    if args.noOfPasses:
        noOfPasses = args.noOfPasses
    print("The number of passes after spline fitting is %s" % noOfPasses) 
    
    noOfBins = 100
    if args.noOfBins:
        noOfBins = args.noOfBins
    print("The number of bins is %s" % noOfBins) 
    
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
    
    visual = False
    if args.visual:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
        visual = True
        print("Graphs will be outputted")
    
    resolution = 0
    if args.resolution:
        resolution = args.resolution
        print("Fixed size option detected... Fast version of FitHiC will be used")
        print("Resolution is %s kb" % (resolution/1000))
    else:
        print("Fixed size data not being used.")  

    interOnly=False
    allReg=False
    if args.chromosome_region:
        if args.chromosome_region == "All":
            allReg=True
            print("All Regions (inter and intrachromosomal)  will be read")
        elif args.chromosome_region == "interOnly":
            interOnly=True
            print("Only interchromosomal interactions will be read")
        elif chromosome_region == "intraOnly":
            interOnly=False
            allReg=False
            print("Only intrachromosmal interactions will be read")
        else:
            print("Invalid Option. Only options are 'All', 'interOnly', or 'intraOnly'") 
    else:
        print("Only intrachromosmal interactions will be read")
     
    print("All arguments processed. Running FitHiC now...")
    print("\n")
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

    baselineIntraChrProb=0  # 1.0/possibleIntraAllCount
    interChrProb=0 # 1.0/possibleInterAllCount

    minObservedGenomicDist=float("inf") 
    maxObservedGenomicDist=0
    maxPossibleGenomicDist=0
    
    #distScaling just avoids overflow - but is necessary for large genomes
    distScaling=10000.0
    
    #RUNBY uses of the following? A: graphs
    toKb=10**-3
    toMb=10**-6
    toProb=10**5

    #intermediate values outputted here
    global logfile
    logfile = os.path.join(outputPath, "fithic.log")
    
    #TODO add visual options including scaling constants

    ##maindic will be generated first using the interactions file only
    mainDic={} # given a distance this dictionary will return [Npairs,TotalContactCount] for only those interactions present in the interactions file
    (mainDic,observedInterAllSum,observedIntraAllSum,observedIntraInRangeSum, biasDic) = read_Interactions(contactCountsFile, biasFile)
    binStats = makeBinsFromInteractions(mainDic, noOfBins, observedIntraInRangeSum)
    
    #Enumerate (fast version) or generate (otherwise) all possible pairs of fragments within the range of interest.
    (mainDic,noOfFrags, maxPossibleGenomicDist, possibleIntraInRangeCount, interChrProb, baselineIntraChrProb)= generate_FragPairs2(binStats, fragsFile, mappThres, resolution, distLowThres, distUpThres)
    
    (binStats,noOfFrags, maxPossibleGenomicDist, possibleIntraInRangeCount, interChrProb, baselineIntraChrProb)= generate_FragPairs3(binStats, fragsFile, bins, resolution)

    

    #read and parse bias values for each locus from ICE or KR normalization output
    if biasFile:
        biasDic = read_biases(biasFile)

    
    #bin the data in desired number of bins, and for each bin, calculate the average genomic distance and average contact probability
    #(x,y,yerr)= calculate_Probabilities(mainDic,resolution,libname+".fithic_pass1", observedIntraInRangeSum,noOfBins, maxPossibleGenomicDist, distScaling)
    (x,y,yerr)= calculate_Probabilities2(mainDic,resolution,libname+".fithic_pass1", observedIntraInRangeSum,noOfBins, maxPossibleGenomicDist, distScaling)
  

    splinefit1st=time.time()
    print("Spline fit Pass 1 starting...")
    #fit a smooth spline to the bin values, and compute and write p values/q values
    splineXinit,splineYinit,splineResidual,isOutlier,splineFDRxinit,splineFDRyinit=fit_Spline(x,y,yerr,options.intersfile,sortedInteractions,biasDic,libname+".spline_pass1",1)
    splinefit1en = time.time()
    print("Spline fit Pass 1 completed. Time took %s" % (splinefit1en-splinefit1st))
    

    if resolution:
        print("Fit-Hi-C completed successfully")

    else:
        ### DO THE NEXT PASSES IF REQUESTED ###
        for i in range(2,2+noOfPasses):
            
            #after removing outliers, rebin the data, and recalculate the average genomic distance and average contact probability
            x,y,yerr=calculate_Probabilities(sortedInteractions,isOutlier,libname+".fithic_pass"+repr(i))
            
            splinefitst=time.time()
            print("Spline fit Pass %s starting..." % i)
            
            #fit a new smooth spline to the recalculated bin values, and recompute and write p values/q values after refinement
            splineX,splineY,splineResidual,isOutlier,splineFDRx,splineFDRy=fit_Spline(x,y,yerr,options.intersfile,sortedInteractions,biasDic,libname+".spline_pass"+repr(i),i)
            splinefiten = time.time()
            print("Spline fit Pass %s completed. Time took %s" % (i, splinefit1en-splinefit1st))



##FUNCTIONS START###

def read_Interactions(contactCountsFile, biasDic):
    mainDic={}
    print("Reading the contact counts file to generate bins...")
    startT = time.time() 

    observedInterAllSum=0 #used
    observedIntraAllSum=0 #used
    observedInterAllCount=0
    observedIntraAllCount=0 #notused
    observedIntraInRangeSum=0 #used
    observedIntraInRangeCount=0 #notused
    minObservedGenomicDist=float("inf")#notused
    maxObservedGenomicDist=0 #notused

    #Loop through every line in the contactCountsFile
    with gzip.open(contactCountsFile, 'r') as f:
        for lines in f:
            ch1,mid1,ch2,mid2,contactCount=lines.split()
           
            #start generating biasDic if not known 
            if biasDic==None:
                if ch1 not in biasDic:
                    biasDic[ch1]={}
                if ch2 not in biasDic:
                    biasDic[ch2]={}
                if mid1 not in biasDic[ch1]:
                    biasDic[ch1][mid1] = 1.0
                if mid2 not in biasDic[ch2]:
                    biasDic[ch2][mid2] = 1.0
        
            #create the interaction
            contactCount=float(contactCount)
            interxn=myUtils.Interaction([ch1, int(mid1), ch2, int(mid2)],distLowThres, distUpThres, contactCount)
            
            if interxn.type=='inter':
                observedInterAllSum += interxn.getCount()
                observedInterAllCount +=1
            else: # any type of intra
                observedIntraAllSum +=interxn.getCount()
                observedIntraAllCount +=1
                if interxn.getType()=='intraInRange':
                    minObservedGenomicDist=min(minObservedGenomicDist,interxn.getDistance())
                    maxObservedGenomicDist=max(maxObservedGenomicDist,interxn.getDistance())
                    if interxn.getDistance() not in mainDic:
                        mainDic[interxn.getDistance()] = [0,0]
                    mainDic[interxn.getDistance()][1]+=interxn.getCount()
                    observedIntraInRangeSum +=interxn.getCount()
                    observedIntraInRangeCount +=1
    endT = time.time()
    print("Interactions file read. Time took %s" % (endT-startT))
    with open(logfile, 'w') as log:
        log.write("Interactions file completed\n")
        log.write("------------------------------------------------------------------------------------\n")
        log.write("Observed, Intra-chr in range: pairs= "+str(observedIntraInRangeCount) +"\t totalCount= "+str(observedIntraInRangeSum)+"\n")
        log.write("Observed, Intra-chr all: pairs= "+str(observedIntraAllCount) +"\t totalCount= "+str(observedIntraAllSum)+"\n")
        log.write("Observed, Inter-chr all: pairs= "+str(observedInterAllCount) +"\t totalCount= "+str(observedInterAllSum)+"\n")
        log.write("Range of observed genomic distances [%d %d]" % (minObservedGenomicDist,maxObservedGenomicDist) + "\n"),
        log.write("\n")
    return (mainDic,observedInterAllSum,observedIntraAllSum,observedIntraInRangeSum, biasDic) # from read_Interactions2

def makeBinsFromInteractions(mainDic,noOfBins, observedIntraInRangeSum):
    with open(logfile, 'a') as log:
        log.write("Making equal occupancy bins\n") 
        log.write("------------------------------------------------------------------------------------\n")
        noPerBin = observedIntraInRangeSum/noOfBins
        log.write("Observed intra-chr read counts in range\t"+repr(observedIntraInRangeSum)+ "\nDesired number of contacts per bin\t" +repr(noPerBin)+",\nNumber of bins\t"+repr(noOfBins)+"\n")
        log.write("\n")

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
    for i in range(len(bins)):
        if i == 0:
            bins[i].append(0)
        else:
            bins[i].append(max(bins[i-1]+1))
    for binIdx in range(len(bins)):
        lb = min(bins[binIdx])
        ub = max(bins[binIdx])
        ##binStats
        #0: range of distances in this bin
        #1: no. of possible pairs w/in this range of distances
        #2: no. of possible pairs INRANGE w/in this range of dists
        #3: sumoverallContactCounts 
        #4: Sumoveralldistances in this bin in millions
        #5: avg CC 
        #6: avg distance
        #7: distances in bin
        binStats[binIdx]=[(lb, ub), 0, 0, 0, 0, 0, 0, bins[binIdx]]
        for dists in bins[binIdx]:
            binStats[binIdx][3]+=mainDic[dists][1]
            binStats[binIdx][4]+=float(dists/distScaling)

    with open(logfile, 'a') as log:
        log.write("Equal occupancy bins generated\n") 
    return binStats

def generate_FragPairs3(binStats, fragsfile, bins, resolution):
    if resolution:
        with open(logfile, 'a') as log:
            log.write("Looping through all possible fragment pairs in-range\n")
            log.write("------------------------------------------------------------------------------------\n"),
    else:
        with open(logfile, 'a') as log:
            log.write("Enumerating all possible fragment pairs in-range\n")
            log.write("------------------------------------------------------------------------------------\n"),
    startT = time.time() 
    
    maxPossibleGenomicDist = 0
    possibleIntraAllCount = 0
    possibleInterAllCount = 0
    possibleIntraInRangeCount = 0
    interChrProb = 0
    baselineIntraChrProb = 0
  
    allFragsDic={}
    with gzip.open(fragsfile,'r') as infile:
        for line in infile:
            words=line.split()
            currChr=words[0]; currMid=words[2];# hitcount=float(words[3]);
            if currChr not in allFragsDic:# and hitcount>=mappThres:
                allFragsDic[currChr]=[]
            allFragsDic[currChr].append(currMid)

    if resolution:
        noOfFrags=0
        maxFrags={}
        
        for ch in allFragsDic:
            maxFrags[ch]=max([int(i)-resolution/2 for i in allFragsDic[ch]])
            noOfFrags+=len(allFragsDic[ch])
            maxPossibleGenomicDist=max(maxPossibleGenomicDist,maxFrags[ch])
       
        for ch in allFragsDic:
            maxFrag=maxFrags[ch]
            n=len(allFragsDic[ch])
            d=0
            binTracker = 0
            possibleIntraInRangeCountPerChr = 0
            for intxnDistance in range(0,maxFrag+1,resolution):
                npairs = n-d
                d+=1
                currBin = binStats[binTracker]
                minOfBin = currBin[0][0]
                maxOfBin = currBin[0][1]
                if minOfBin<=intxnDistance<=maxOfBin:
                    pass
                else:
                    binTracker+=1
                    if binTracker not in binStats:
                        binTracker-=1
                currBin[1]+=npairs
                currBin[4]+=float(intxnDistance/distScaling)
                if myUtils.in_range_check(intxnDistance,distLowThres,distUpThres):
                    currBin[2]+=npairs
                    possibleIntraInRangeCountPerChr += 1
            possibleInterAllCount+=n*(noOfFrags-n)
            possibleIntraAllCount+=(n*(n+1))/2 # n(n-1) if excluding self
            with open(logfile, 'a') as log:
                log.write("Chromosome " +repr(ch) +",\t"+str(templen) +" mappable fragments, \t"+str(possibleIntraInRangeCountPerChr)\
                +" possible intra-chr fragment pairs in range,\t" + str((noOfFrags-templen)*templen) +" possible inter-chr fragment pairs\n")
            posibleIntraInRangeCount += possibleIntraInRangeCountPerChr
        possibleInterAllCount/=2
        interChrProb=1.0/possibleInterAllCount
        baselineIntraChrProb=1.0/possibleIntraAllCount
        
    else:
        noOfFrags = 0
        for ch in allFragsDic:
            noOfFrags += len(allFragsDic[ch])

        for ch in allFragsDic: 
            countIntraPairs = 0
            fragsPerChr = sorted(allFragsDic[ch])
            templen = len(fragsPerChr)
            possibleInterAllCount += (noOfFrags-templen)*templen
            binTracker = 0
            possibleIntraInRangeCountPerChr = 0
            for x in range(templen):
                d = 0
                for y in range(x+1,templen):
                    intxnDistance = abs(float(fragsPerChr[x])-float(fragsPerChr[y]))
                    maxPossibleGenomicDist = max(maxPossibleGenomicDist, intxnDistance)
                    npairs = templen-d
                    d+=1
                    currBin = binStats[binTracker]
                    minOfBin = currBin[0][0]
                    maxOfBin = currBin[0][1]
                    if minOfBin<=intxnDistance<=maxOfBin:
                        pass
                    else:
                        binTracker+=1
                        if binTracker not in binStats:
                            binTracker-=1
                    currBin[1]+=npairs
                    currBin[4]+=float(intxnDistance/distScaling)
                    if myUtils.in_range_check(intxnDistance,distLowThres,distUpThres):
                        currBin[2]+=npairs
                        possibleIntraInRangeCountPerChr += 1
                    possibleIntraAllCount += 1
            with open(logfile, 'a') as log:
                log.write("Chromosome " +repr(ch) +",\t"+str(templen) +" mappable fragments, \t"+str(possibleIntraInRangeCountPerChr)\
                +" possible intra-chr fragment pairs in range,\t" + str((noOfFrags-templen)*templen) +" possible inter-chr fragment pairs\n")
            possibleIntraInRangeCount += possibleIntraInRangeCountPerChr
        possibleInterAllCount/=2
        interChrProb=1.0/possibleInterAllCount
        baselineIntraChrProb=1.0/possibleIntraAllCount
    endT = time.time()
    print("Fragments file read. Time took %s" % (endT-startT))

    with open(logfile, 'a') as log:
        log.write("Number of all fragments= %s\n" % (noOfFrags))
        log.write("Possible, Intra-chr in range: pairs= %s \n" % (possibleIntraInRangeCount))
        log.write("Possible, Intra-chr all: pairs= %s \n" % (possibleIntraAllCount)) 
        log.write("Possible, Inter-chr all: pairs= %s \n" % (possibleInterAllCount)) log.write("Desired genomic distance range   [%d %s] \n" % (distLowThres,distUpThres)),
        log.write("Range of possible genomic distances  [0  %d] \n" % (maxPossibleGenomicDist)),

    return (binStats,noOfFrags, maxPossibleGenomicDist, possibleIntraInRangeCount, interChrProb, baselineIntraChrProb) # return from generate_FragPairs3


def read_biases(infilename):
    print("\n\nReading biases. \n")
    biasDic={}
    
    rawBiases=[]
    infile =gzip.open(infilename, 'r')
    for line in infile:
        words=line.rstrip().split()
        chr=words[0]; midPoint=int(words[1]); bias=float(words[2])
        
        #NEW remove this
        #if bias!=1.0:
        #   rawBiases.append(bias)
    infile.close()
    #sys.stderr.write("\n\nReading ICE biases. \n")
    botQ,med,topQ=mquantiles(rawBiases,prob=[0.05,0.5,0.95])
    with open(logfile, 'a') as log:
        log.write("5th quantile of biases: "+str(botQ)+"\n")
        log.write("50th quantile of biases: "+str(med)+"\n")
        log.write("95th quantile of biases: "+str(topQ)+"\n")

    #m,v=myStats.meanAndVariance(rawBiases)
    #sd=math.sqrt(v)
    #sys.stderr.write(str(m)+"\t"+str(v)+"\t"+str(sd)+"\n")
    
    #normFactor=sum(rawBiases)/len(rawBiases)
    infile =gzip.open(infilename, 'r')
    totalC=0
    discardC=0
    for line in infile: words=line.rstrip().split()
        chr=words[0]; midPoint=int(words[1]); bias=float(words[2])
        if bias<0.5:
            bias=-1 #botQ
            discardC+=1
        elif bias>2:
            bias=-1 #topQ
            discardC+=1
        totalC+=1
        if chr not in biasDic:
            biasDic[chr]={}
        if midPoint not in biasDic[chr]:
            biasDic[chr][midPoint]=bias
    infile.close()
    with open(logfile, 'a') as log:
        log.write("Out of " + str(totalC) + " loci " +str(discardC) +" were discarded with biases not in range [0.5 2]\n\n" )
    print("Bias file read")
    return biasDic # from read_biases

def calculateProbabilities2(binStats,resolution,outfilename,observedIntraInRangeSum,noOfBins, maxPossibleGenomicDist, distScaling):
    with open(logfile, 'a') as log:
        log.write("\nCalculating probability means and standard deviations of contact counts\n"),
        log.write("------------------------------------------------------------------------------------\n"),
        log.write("observed intra-chr read counts in range\t"+repr(observedIntraInRangeSum)+ ",\tdesired number of contacts per bin\t" +repr(desiredPerBin)+",\tnumber of bins\t"+repr(noOfBins)+"\n"),
    
    if resolution:
        outfile=open(outfilename+'.res'+str(resolution)+'.txt', 'w') 
    else:
        outfile=open(outfilename+'.txt', 'w')
   
    x = []
    y = []
    yerr = []
    pairCounts=[]
    interactionTotals=[]

    for i in range(len(binStats)):
        currBin = binStats[i]
        sumCC = currBin[3]
        sumDistB4Scaling = currBin[4]
        possPairsInRange = currBin[2] #RUNBY: currBin[1]?
        avgCC = (1.0*sumCC/possPairsInRange)/observedIntraInRangeSum
        avgDist = distScaling*(sumDistB4Scaling/possPairsInRange)
        currBin[5]=avgCC
        currBin[6]=avgDist
        y.append(avgCC)
        x.append(avgDist)
        yerr.append(0) #TODO FIX
        pairCounts.append(possPairsInRange)
        interactionTotals.append(sumCC)
        
    print("Writing %s" % outfile"\n")
    outfile.write("avgGenomicDist\tcontactProbability\tstandardError\tnoOfLocusPairs\ttotalOfContactCounts\n")
    for i in range(len(x)):
        outfile.write("%d" % x[i] + "\t"+"%.2e" % y[i]+ "\t" + "%.2e" % yerr[i] + "\t" +"%d" % pairCounts[i] + "\t" +"%d" % interactionTotals[i]+"\n")
    outfile.close()
    return [x,y,yerr] # from calculateProbabilities2


def fit_Spline(mainDic,x,y,yerr,infilename,outfilename,biasDic):
    with open(logfile, 'a') as log:
        log.write("\nFit a univariate spline to the probability means\n"),
        log.write("------------------------------------------------------------------------------------\n"),
    #print("baseline intra-chr probability: " + repr(baselineIntraChrProb)+ "\n"),

    # maximum residual allowed for spline is set to min(y)^2
    splineError=min(y)*min(y) 

    # use fitpack2 method -fit on the real x and y from equal occupancy binning
    ius = UnivariateSpline(x, y, s=splineError)

    #### POST-PROCESS THE SPLINE TO MAKE SURE IT'S NON-INCREASING
    ### NOW I DO THIS BY CALLING AN R function CALLED MONOREG
    ### This does the isotonic regression using option antitonic to make sure
    ### I get monotonically decreasing probabilites with increasion genomic distance

    tempMaxX=max(x)
    tempMinX=min(x)
    tempList=sorted([dis for dis in mainDic])
    splineX=[]
    ### The below for loop will make sure nothing is out of range of [min(x) max(x)]
    ### Therefore everything will be within the range where the spline is defined
    for i in tempList:
        if tempMinX<=i and i<=tempMaxX:
            splineX.append(i)
    # END for
    splineY=ius(splineX)

    # R vector format
    rSplineX=ro.FloatVector(splineX)
    rSplineY=ro.FloatVector(splineY)
    rMonoReg=ro.r['monoreg']
    # do the antitonic regression
    allRres=rMonoReg(rSplineX,rSplineY,type="antitonic")
    rNewSplineY=allRres[3]
    # convert data back to Python format
    newSplineY=[]
    diff=[]
    diffX=[]
    for i in range(len(rNewSplineY)):
        newSplineY.append(rNewSplineY[i])
        if (splineY[i]-newSplineY[i]) > 0:
            diff.append(splineY[i]-newSplineY[i])
            diffX.append(splineX[i])
        # END if
    # END for

    ### Now newSplineY holds the monotonic contact probabilities
    residual =sum([i*i for i in (y - ius(x))])

    ### Now plot the results
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(2,1,1)
    plt.title('Univariate spline fit to the output of equal occupancy binning. \n Residual= %e' % (residual),size='small')
    plt.plot([i/1000.0 for i in x], [i*100000 for i in y], 'ro', label="Means")
    #plt.plot([i/1000.0 for i in xi], [i*100000 for i in yi],'g-',label="Spline fit")
    plt.plot([i/1000.0 for i in splineX], [i*100000 for i in newSplineY],'g-',label="Spline fit")
    #plt.plot([i/1000.0 for i in x], [normalizedInterChrProb*100000 for i in x],'k-',label="Random intra-chromosomal")
    #plt.plot([i/1000.0 for i in x], [interChrProb*100000 for i in x],'b-',label="Inter-chromosomal")
    plt.ylabel('Probability (1e-5)')
    plt.xlabel('Genomic distance (kb)')
    plt.xlim([min(x)/1000.0,max(x)/1000.0])
    ax.legend(loc="upper right")

    ax = fig.add_subplot(2,1,2)
    plt.loglog(splineX,newSplineY,'g-')
    #plt.loglog(xi, yi, 'g-') 
    plt.loglog(x, y, 'r.')  # Data
    #plt.loglog(x, [normalizedInterChrProb for i in x],'k-')
    #plt.loglog(x, [interChrProb for i in x],'b-')
    plt.ylabel('Probability (log scale)')
    plt.xlabel('Genomic distance (log scale)')
    #plt.xlim([20000,100000])
    plt.xlim([min(x),max(x)])
    plt.savefig(outfilename+'.res'+str(resolution)+'.png')
    sys.stderr.write("Plotting %s" % outfilename + ".png\n")

    # NOW write the calculated pvalues and corrected pvalues in a file
    infile =gzip.open(infilename, 'r')
    intraInRangeCount=0
    intraOutOfRangeCount=0
    intraVeryProximalCount=0
    interCount=0
    discardCount=0
    print("lower bound on mid-range distances  "+ repr(distLowThres) + ", upper bound on mid-range distances  " + repr(distUpThres) +"\n"),
    p_vals=[]
    q_vals=[]
    for line in infile:
        words=line.rstrip().split()
        interxn=myUtils.Interaction([words[0], int(words[1]), words[2], int(words[3])])
        interxn.setCount(float(words[4]))
        chr1=words[0]
        chr2=words[2]
        midPoint1=int(words[1])
        midPoint2=int(words[3])
        
        bias1=1.0; bias2=1.0;  # assumes there is no bias to begin with
        # if the biasDic is not null sets the real bias values
        if len(biasDic)>0:
            if chr1 in biasDic and midPoint1 in biasDic[chr1]:
                bias1=biasDic[chr1][midPoint1]
            if chr2 in biasDic and midPoint2 in biasDic[chr2]:
                bias2=biasDic[chr2][midPoint2]
    
        if bias1==-1 or bias2==-1:
            p_val=1.0
            discardCount+=1
        elif interxn.type=='intra':
            if interxn.getType(distLowThres,distUpThres)=='intraInRange':
                # make sure the interaction distance is covered by the probability bins
                distToLookUp=max(interxn.distance,min(x))
                distToLookUp=min(distToLookUp,max(x))
                i=min(bisect.bisect_left(splineX, distToLookUp),len(splineX)-1)
                prior_p=newSplineY[i]*(bias1*bias2) # biases added in the picture
                p_val=scsp.bdtrc(interxn.hitCount-1,observedIntraInRangeSum,prior_p)
                intraInRangeCount +=1
            elif interxn.getType(distLowThres,distUpThres)=='intraShort':
                prior_p=1.0
                p_val=1.0
                intraVeryProximalCount +=1
            elif interxn.getType(distLowThres,distUpThres)=='intraLong':
                ## out of range distance
                ## use the prior of the baseline intra-chr interaction probability
                prior_p=baselineIntraChrProb*(bias1*bias2)  # biases added in the picture
                p_val=scsp.bdtrc(interxn.hitCount-1,observedIntraAllSum,prior_p)
                intraOutOfRangeCount +=1
            # END if
        else: # inter
            #prior_p=normalizedInterChrProb
            prior_p=interChrProb*(bias1*bias2) # biases added in the picture
            ############# THIS HAS TO BE interactionCount-1 ##################
            p_val=scsp.bdtrc(interxn.hitCount-1,observedInterAllSum,prior_p)
            interCount +=1
        #
        p_vals.append(p_val)

    # END for
    infile.close()

    # Do the BH FDR correction
    ## This below line was INCORRECTLY overcorrecting pvalues - commented out June 13th, 2017 by Ferhat
    #q_vals=myStats.benjamini_hochberg_correction(p_vals, possibleInterAllCount+possibleIntraAllCount)
    q_vals=myStats.benjamini_hochberg_correction(p_vals, possibleIntraInRangeCount)
    #print("possibleIntraInRangeCount " + repr(possibleIntraInRangeCount)+"\n"),

    infile =gzip.open(infilename, 'r')
    outfile =gzip.open(outfilename+'.res'+str(resolution)+'.significances.txt.gz', 'w')
    print("Writing p-values and q-values to file %s" % outfilename + ".significances.txt\n"),
    print("Number of pairs discarded due to bias not in range [0.5 2]\n"),
    outfile.write("chr1\tfragmentMid1\tchr2\tfragmentMid2\tcontactCount\tp-value\tq-value\tbias1\tbias2\n")
    count=0
    for line in infile:
        words=line.rstrip().split()
        chr1=words[0]
        midPoint1=int(words[1])
        chr2=words[2]
        midPoint2=int(words[3])
        interactionCount=int(words[4])
        p_val=p_vals[count]
        q_val=q_vals[count]
        bias1=biasDic[chr1][midPoint1]
        bias2=biasDic[chr2][midPoint2]

        #if chr1==chr2: # intra
        #   interactionDistance=abs(midPoint1-midPoint2) # dist
        #   if myUtils.in_range_check(interactionDistance,distLowThres,distUpThres):
        #       outfile.write("%s\t%d\t%s\t%d\t%d\t%e\t%e\n" % (str(chr1),midPoint1,str(chr2),midPoint2,interactionCount,p_val,q_val))
        #else:
        #   outfile.write("%s\t%d\t%s\t%d\t%d\t%e\t%e\n" % (str(chr1),midPoint1,str(chr2),midPoint2,interactionCount,p_val,q_val))

        outfile.write("%s\t%d\t%s\t%d\t%d\t%e\t%e\t%.3f\t%.3f\n" % (str(chr1),midPoint1,str(chr2),midPoint2,interactionCount,p_val,q_val,bias1,bias2))
        count+=1
    # END for - printing pvals and qvals for all the interactions
    outfile.close()
    infile.close()
    return [splineX, newSplineY, residual] # from fit_Spline

def calculate_Probabilities(binStats,resolution,outfilename,observedIntraInRangeSum,noOfBins, maxPossibleGenomicDist, distScaling): #TODO fix resolution
    desiredPerBin=(observedIntraInRangeSum)/noOfBins
    with open(logfile, 'a') as log:
        log.write("\nCalculating probability means and standard deviations by equal-occupancy binning of contact counts\n"),
        log.write("------------------------------------------------------------------------------------\n"),
        log.write("observed intra-chr read counts in range\t"+repr(observedIntraInRangeSum)+ ",\tdesired number of contacts per bin\t" +repr(desiredPerBin)+",\tnumber of bins\t"+repr(noOfBins)+"\n"),
    
    if resolution:
        outfile=open(outfilename+'.res'+str(resolution)+'.txt', 'w') 
    else:
        outfile=open(outfilename+'.txt', 'w')
    
    # the following five lists will be the print outputs
    x=[] # avg genomic distances of bins
    y=[] # avg interaction probabilities of bins
    yerr=[] # stderrs of bins
    pairCounts=[] # number of pairs in bins
    interactionTotals=[] # number of interactions (reads) in bins
    interactionTotalForBinTermination=0
    n=0 # bin counter so far
    totalInteractionCountSoFar=0
    distsToGoInAbin=[]
    binFull=0
    for i in range(0,maxPossibleGenomicDist+1,resolution):
        totalInteractionCountSoFar+=mainDic[i][1]
        if myUtils.in_range_check(i,distLowThres,distUpThres)==False:
            continue
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
        # if adding the next bin will fill the bin
        else:
            distsToGoInAbin.append(i)
            interactionTotalForBinTermination+=mainDic[i][1]
        #
        if binFull==1:
            noOfPairsForBin=0
            interactionTotalForBin=0
            avgDistance=0
            # dynamically update the desiredPerBin after each bin is full
            n+=1
            if n<noOfBins:
                desiredPerBin=1.0*(observedIntraInRangeSum-totalInteractionCountSoFar)/(noOfBins-n)
            se_p=0 # for now I'm not worrying about error etc.
            for b in distsToGoInAbin:
                noOfPairsForBin+=mainDic[b][0]
                interactionTotalForBin+=mainDic[b][1]
                avgDistance+=1.0*mainDic[b][0]*(b/distScaling)
            #
            meanProbabilityObsv=(1.0*interactionTotalForBin/noOfPairsForBin)/observedIntraInRangeSum
            avgDistance=distScaling*(avgDistance/noOfPairsForBin)
            # append this bin
            x.append(float(avgDistance))
            y.append(float(meanProbabilityObsv))
            yerr.append(float(se_p))
            pairCounts.append(noOfPairsForBin)
            interactionTotals.append(interactionTotalForBin)
            
            print "%d" % n+ "\t" + "%f" % avgDistance + "\t"+"%.2e" % meanProbabilityObsv + "\t"\
                + "%.2e" % se_p +"\t" +"%d" % noOfPairsForBin +"\t" +"%d" % interactionTotalForBin
            # reset counts
            interactionTotalForBinTermination=0
            binFull=0
            distsToGoInAbin=[]
        # END if
    # END for
    with open(logfile, 'a') as log:

        log.write("Writing equal-occupancy binning results to %s" % outfilename + ".txt\n"),
    outfile.write("avgGenomicDist\tcontactProbability\tstandardError\tnoOfLocusPairs\ttotalOfContactCounts\n")
    for i in range(len(x)):
        outfile.write("%d" % x[i] + "\t"+"%.2e" % y[i]+ "\t" + "%.2e" % yerr[i] + "\t" +"%d" % pairCounts[i] + "\t" +"%d" % interactionTotals[i]+"\n")
    outfile.close()
    return [x,y,yerr] # from calculate_Probabilities




if __name__ == "__main__":
    main()
