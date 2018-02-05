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
                        #implement this TODO


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
                      distance range (unit: base pairs). DEFAULT no limit. 
                      STRONGLY suggested to have a limit for large genomes,
                      such as human/mouse. ex. '1000000, 5000000, etc.'",
                      required=False)

    parser.add_argument("-L", "--lowerbound", dest="distLowThres", type=int,
                      help="OPTIONAL: lower bound on the intra-chromosomal \
                      distance range (unit: base pairs). DEFAULT no limit.
                      Suggested limit is 2x the resolution of the input files",
                      required=False)

    parser.add_argument("-v", "--visual", action="store_true", dest="visual",\
                      help="OPTIONAL: use this flag for generating plots. \
                      DEFAULT is False.", required=False)
    
    parser.add_argument("-V", "--version", action="version",version=versionStr)

    parser.add_argument("-x", "--chromosome_region", dest="chromosome_region", 
                      help="OPTIONAL: use this flag to determine which chromosomal \
                      regions to study (intraOnly, interOnly, All) \
                      DEFAULT is intraOnly", required=False)
    
    parser.add_argument("-r", "--resolution", dest="resolution", type=int,  
                      help="OPTIONAL: Use this option if the fragments file \
                      is fixed size binned data). DEFAULT is None", required=False)
    return parser.parse_args()

def main():
    args = parse_args(sys.argv[1:])
    

    #TODO add '-x' parsing
    #TODO sort by coded
    #TODO add error/warning messages
    fragsFile = args.fragsfile
    contactCountsFile = args.intersfile
    outputPath = args.outdir
    
    noOfPasses = 0
    if args.noOfPasses:
        noOfPasses = args.noOfPasses
    
    libName = "FitHiC"
    if args.libname:
        libName = args.libname
    
    noOfBins = 100
    if args.noOfBins:
        noOfBins = args.noOfBins
    
    biasFile = None
    if args.biasfile:
        biasFile = args.biasfile

    mappThres = 1
    if args.mappabilityThreshold:
        mappThres = args.mappabilityThreshold

    distUpThres = 0 #TODO change to inf.
    distLowThres = 0
    if args.distUpThres:
        distUpThres = args.distUpThres
    if args.distLowThres:
        distLowThres = args.distLowThres

    visual = False
    if args.visual:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
        visual = True

    resolution = 0
    if args.resolution:
        resolution = args.resolution
    #TODO add summary print
    
    #TODO end
    #All arguments processed
   

    possibleIntraInRangeCount=0 # count of all possible inter-chr fragment pairs
    observedIntraInRangeCount=0
    observedIntraInRangeSum=0
    # intra-chromosomal contacts
    possibleIntraAllCount=0 # count of all possible intra-chr fragment pairs
    observedIntraAllCount=0
    observedIntraAllSum=0
    # inter-chromosomal contacts
    possibleInterAllCount=0 # count of all possible inter-chr fragment pairs
    observedInterAllCount=0
    observedInterAllSum=0

    baselineIntraChrProb=0	# 1.0/possibleIntraAllCount
    interChrProb=0 #  1.0/possibleInterAllCount

    minObservedGenomicDist=500000000 # some number bigger than the biggest chromosome length #TODO change to inf.
    maxObservedGenomicDist=0
    maxPossibleGenomicDist=0
    
    #distScaling just avoids overflow - but is necessary for large genomes
    distScaling=10000.0
    
    #RUNBY uses of the following?
    toKb=10**-3
    toMb=10**-6
    toProb=10**5
    #TODO add visual options including scaling constants


	mainDic={} # given a distance this dictionary will return [Npairs,TotalContactCount] including all possible pairs (even if not in input file)

	#Enumerate (fast version) or generate (otherwise) all possible pairs of fragments within the range of interest.
    (mainDic,noOfFrags, maxPossibleGenomicDist, possibleIntraInRangeCount, interChrProb, baselineIntraChrProb)= generate_FragPairs(fragsFile, resolution))
    
    #read and parse bias values for each locus from ICE or KR normalization output
    if biasFile:
        biasDic = read_biases(biasFile)
    
    #reading all interaction/contact counts from fithic formatted input file
    (mainDic,observedInterAllSum,observedIntraAllSum,observedIntraInRangeSum, biasDic)=read_All_Interactions(mainDic,contactCountsFile,noOfFrags)
   

    #TODO probably will not need this 
    print("\n\t\tSPLINE FIT PASS 1 of %s \n" % (noOfPasses+1)))
    
    #bin the data in desired number of bins, and for each bin, calculate the average genomic distance and average contact probability
    (x,y,yer)= calculate_Probabilities(mainDic,resolution,outfilename, observedIntraInRangeSum,noOfBins, maxPossibleGenomicDist, distScaling)
   
    #fit a smooth spline to the bin values, and compute and write p values/q values
    splineXinit,splineYinit,splineResidual,isOutlier,splineFDRxinit,splineFDRyinit=fit_Spline(x,y,yerr,options.intersfile,sortedInteractions,biasDic,libname+".spline_pass1",1)


    #TODO if statement if not fast version and noofpasses is nonzero
    ### DO THE NEXT PASSES IF REQUESTED ###
    for i in range(2,2+noOfPasses):
        
        #after removing outliers, rebin the data, and recalculate the average genomic distance and average contact probability
        x,y,yerr=calculate_Probabilities(sortedInteractions,isOutlier,libname+".fithic_pass"+repr(i))
        
        
        #fit a new smooth spline to the recalculated bin values, and recompute and write p values/q values after refinement
        splineX,splineY,splineResidual,isOutlier,splineFDRx,splineFDRy=fit_Spline(x,y,yerr,options.intersfile,sortedInteractions,biasDic,libname+".spline_pass"+repr(i),i)



##FUNCTIONS START###



def generate_FragPairs(fragsfile, mappThres, resolution):
	print("\nEnumerating all possible fragment pairs in-range\n"),
	print("------------------------------------------------------------------------------------\n"),
	mainDic={} 
    maxPossibleGenomicDist = 0
    possibleIntraAllCount = 0
	possibleInterAllCount = 0
	possibleIntraInRangeCount = 0
	interChrProb = 0
	baselineIntraChrProb = 0
	#badFrags=[]
	allFragsDic={}
	#allFragsDicReverse={}
	infile=gzip.open(fragsfile,'r')
	indx=0
	for line in infile:
		words=line.split()
		currChr=words[0]; currMid=words[2]; mapp=float(words[4]); #RUNBY fileformat/mapp
		if currChr not in allFragsDic and mappThres<=mapp:
			allFragsDic[currChr]={}
		allFragsDic[currChr][currMid]=indx
	#	allFragsDicReverse[indx]=[currChr,currMid]
		#if mapp<=lowMappThres:
		#	badFrags.append(indx)
		indx+=1
	infile.close()

	noOfFrags=0
	maxFrags={}
	for ch in allFragsDic:
		maxFrags[ch]=max([int(i)-resolution/2 for i in allFragsDic[ch]])
		noOfFrags+=len(allFragsDic[ch])
		maxPossibleGenomicDist=max(maxPossibleGenomicDist,maxFrags[ch])
	#print badFrags

	for i in range(0,maxPossibleGenomicDist+1,resolution):
		mainDic[i]=[0,0]
    #TODO 'set' stage for resolution, first read all interactions and then create bins with respect to those with nonzero interactions. After the fact redo enumeration. make bins first and then sort by each frag. Enumerateion is already sorting the data. write function that first reads interaction counts from inter counts file, do binning on that, and then loop through fragfile. If fast do tiny loop, if not fast then do iterations. data function that sees if 0 contact. 
    #already finished binning withrespect to those dist., for each dist. how many counts and then look at the bin, how many dist. does this bin span, those 3 you sum up and put there

#readprobabilities compute bin contact counts, now call function that behaves diff. whether fast/not to get denominators and then do division right after

	for ch in allFragsDic:
		maxFrag=maxFrags[ch]
		n=len(allFragsDic[ch])
		d=0
		for i in range(0,maxFrag+1,resolution):
			mainDic[i][0]+=n-d
			d+=1
		#
		possibleInterAllCount+=n*(noOfFrags-n)
		possibleIntraAllCount+=(n*(n+1))/2 # n(n-1) if excluding self
	#
	possibleInterAllCount/=2
	interChrProb=1.0/possibleInterAllCount
	baselineIntraChrProb=1.0/possibleIntraAllCount
	
	for i in range(0,maxPossibleGenomicDist+1,resolution):
		if myUtils.in_range_check(i,distLowThres,distUpThres):
			possibleIntraInRangeCount+=mainDic[i][0]
		#print str(i)+"\t"+str(mainDic[i][0])

	print("Number of all fragments= "+str(noOfFrags)+"\t resolution= "+ str(resolution))
	print("Possible, Intra-chr in range: pairs= "+str(possibleIntraInRangeCount))
	print("Possible, Intra-chr all: pairs= "+str(possibleIntraAllCount)) 
	print("Possible, Inter-chr all: pairs= "+str(possibleInterAllCount))
	print("Desired genomic distance range	[%d %d]" % (distLowThres,distUpThres) + "\n"),
	print("Range of possible genomic distances	[0	%d]" % (maxPossibleGenomicDist) + "\n"),
	
    maxPossibleGenomicDist = 0
	possibleIntraInRangeCount = 0
	interChrProb = 0
	baselineIntraChrProb = 0
	#badFrags=[]
	allFragsDic={}

	return (mainDic,noOfFrags, maxPossibleGenomicDist, possibleIntraInRangeCount, interChrProb, baselineIntraChrProb) # return from generate_FragPairs

def read_biases(infilename):
	sys.stderr.write("\n\nReading ICE biases. \n")
	biasDic={}
	
	rawBiases=[]
	infile =gzip.open(infilename, 'r')
	for line in infile:
		words=line.rstrip().split()
		chr=words[0]; midPoint=int(words[1]); bias=float(words[2])
		
        #NEW remove this
        #if bias!=1.0:
		#	rawBiases.append(bias)
	infile.close()
	#sys.stderr.write("\n\nReading ICE biases. \n")
	botQ,med,topQ=mquantiles(rawBiases,prob=[0.05,0.5,0.95])
	sys.stderr.write("5th quantile of biases: "+str(botQ)+"\n")
	sys.stderr.write("50th quantile of biases: "+str(med)+"\n")
	sys.stderr.write("95th quantile of biases: "+str(topQ)+"\n")

	#m,v=myStats.meanAndVariance(rawBiases)
	#sd=math.sqrt(v)
	#sys.stderr.write(str(m)+"\t"+str(v)+"\t"+str(sd)+"\n")
	
	#normFactor=sum(rawBiases)/len(rawBiases)
	infile =gzip.open(infilename, 'r')
	totalC=0
	discardC=0
	for line in infile:
		words=line.rstrip().split()
		chr=words[0]; midPoint=int(words[1]); bias=float(words[2])
		# extra conditions
		#if bias<(botQ/2.0):
		if bias<0.5:
			bias=-1 #botQ
			discardC+=1
		elif bias>2:
			bias=-1 #topQ
			discardC+=1
		#
		totalC+=1
		if chr not in biasDic:
			biasDic[chr]={}
		if midPoint not in biasDic[chr]:
			biasDic[chr][midPoint]=bias
	infile.close()
	sys.stderr.write("Out of " + str(totalC) + " loci " +str(discardC) +" were discarded with biases not in range [0.5 2]\n\n" )

	return biasDic # from read_biases

def	read_All_Interactions(mainDic,contactCountsFile,noOfFrags,biasDic):
	print("\nReading all the contact counts\n"),
	print("------------------------------------------------------------------------------------\n"),

	observedInterAllSum=0 #used
	observedInterAllCount=0 #notused
	observedIntraAllSum=0 #used
	observedIntraAllCount=0 #notused
	observedIntraInRangeSum=0 #used
	observedIntraInRangeCount=0 #notused
	minObservedGenomicDist=0 #notused
	maxObservedGenomicDist=0 #notused

	#Xvals=[]
	#Xindices=[]
	#for i in range(noOfFrags):
	#	Xvals.append([])
	#	Xindices.append([])
	##
	infile=gzip.open(contactCountsFile,'r')
	count=0
	for line in infile:
		ch1,mid1,ch2,mid2,contactCount=line.split()

        if biasDic==None:
            if chr1 not in biasDic:
                biasDic[chr1]={}
            if chr2 not in biasDic:
                biasDic[chr2]={}
            if mid1 not in biasDic[chr1]:
                biasDic[chr1][mid1] = 1.0
            if mid2 not in biasDic[chr2]:
                biasDic[chr2][mid2] = 1.0
		### FIXME: this part will need to be fixed for human etc
		#ch1='chr'+ch1
		#ch2='chr'+ch2
		contactCount=float(contactCount)
		interxn=myUtils.Interaction([ch1, int(mid1), ch2, int(mid2)])
		interxn.setCount(contactCount)
		count+=1

		if count%1000000==0:
			print count
		if interxn.type=='inter':
			observedInterAllSum +=interxn.hitCount
			observedInterAllCount +=1
		else: # any type of intra
			observedIntraAllSum +=interxn.hitCount
			observedIntraAllCount +=1
			if interxn.getType(distLowThres,distUpThres)=='intraInRange':
				minObservedGenomicDist=min(minObservedGenomicDist,interxn.distance)
				maxObservedGenomicDist=max(maxObservedGenomicDist,interxn.distance)
				if interxn.distance in mainDic:
					mainDic[interxn.distance][1]+=contactCount
				observedIntraInRangeSum +=interxn.hitCount
				observedIntraInRangeCount +=1
		# END else
	#	indx1=allFragsDic[ch1][mid1]
	#	indx2=allFragsDic[ch2][mid2]
		#print str(indx1)+"\t"+str(indx2)
	#	Xvals[indx1].append(contactCount)
	#	Xindices[indx1].append(indx2)
	#	Xvals[indx2].append(contactCount)
	#	Xindices[indx2].append(indx1)
	# END for
	infile.close()
	print("Observed, Intra-chr in range: pairs= "+str(observedIntraInRangeCount) +"\t totalCount= "+str(observedIntraInRangeSum))
	print("Observed, Intra-chr all: pairs= "+str(observedIntraAllCount) +"\t totalCount= "+str(observedIntraAllSum))
	print("Observed, Inter-chr all: pairs= "+str(observedInterAllCount) +"\t totalCount= "+str(observedInterAllSum))
	print("Range of observed genomic distances [%d %d]" % (minObservedGenomicDist,maxObservedGenomicDist) + "\n"),

	return (mainDic,observedInterAllSum,observedIntraAllSum,observedIntraInRangeSum, biasDic) # from read_All_Interactions

def calculate_Probabilities(mainDic,resolution,outfilename, observedIntraInRangeSum,noOfBins, maxPossibleGenomicDist, distScaling): #TODO fix resolution
	print("\nCalculating probability means and standard deviations by equal-occupancy binning of contact counts\n"),
	print("------------------------------------------------------------------------------------\n"),
	outfile=open(outfilename+'.res'+str(resolution)+'.txt', 'w') #TODO change this to open/with

	## total interaction count to put on top of the plot
	#totalInteractionCountForPlot=0
	#for i in range(0,maxPossibleGenomicDist+1,resolution):
	#	totalInteractionCountForPlot += mainDic[i][1]
	#totalInteractionCountForPlot/=2

	desiredPerBin=(observedIntraInRangeSum)/noOfBins
	print("observed intra-chr read counts in range\t"+repr(observedIntraInRangeSum)+ ",\tdesired number of contacts per bin\t" +repr(desiredPerBin)+",\tnumber of bins\t"+repr(noOfBins)+"\n"),

	# the following five lists will be the print outputs
	x=[] # avg genomic distances of bins
	y=[] # avg interaction probabilities of bins
	yerr=[] # stderrs of bins
	pairCounts=[] # number of pairs in bins
	interactionTotals=[] # number of interactions (reads) in bins
	interactionTotalForBinTermination=0
	n=0 # bin counter so far
	totalInteractionCountSoFar=0
	#observedIntraInRangeSum
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
	print("Writing equal-occupancy binning results to %s" % outfilename + ".txt\n"),
	outfile.write("avgGenomicDist\tcontactProbability\tstandardError\tnoOfLocusPairs\ttotalOfContactCounts\n")
	for i in range(len(x)):
		outfile.write("%d" % x[i] + "\t"+"%.2e" % y[i]+ "\t" + "%.2e" % yerr[i] + "\t" +"%d" % pairCounts[i] + "\t" +"%d" % interactionTotals[i]+"\n")
	outfile.close()
	return [x,y,yerr] # from calculate_Probabilities


def fit_Spline(mainDic,x,y,yerr,infilename,outfilename,biasDic):
	print("\nFit a univariate spline to the probability means\n"),
	print("------------------------------------------------------------------------------------\n"),
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
		#	interactionDistance=abs(midPoint1-midPoint2) # dist
		#	if myUtils.in_range_check(interactionDistance,distLowThres,distUpThres):
		#		outfile.write("%s\t%d\t%s\t%d\t%d\t%e\t%e\n" % (str(chr1),midPoint1,str(chr2),midPoint2,interactionCount,p_val,q_val))
		#else:
		#	outfile.write("%s\t%d\t%s\t%d\t%d\t%e\t%e\n" % (str(chr1),midPoint1,str(chr2),midPoint2,interactionCount,p_val,q_val))

		outfile.write("%s\t%d\t%s\t%d\t%d\t%e\t%e\t%.3f\t%.3f\n" % (str(chr1),midPoint1,str(chr2),midPoint2,interactionCount,p_val,q_val,bias1,bias2))
		count+=1
	# END for - printing pvals and qvals for all the interactions
	outfile.close()
	infile.close()
	return [splineX, newSplineY, residual] # from fit_Spline


if __name__ == "__main__":
    main()
