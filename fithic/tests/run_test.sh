#!/bin/bash -ex

set -o nounset
set -o pipefail
set -o errexit 

BINDIR=../bin
DATADIR=data
inI=$DATADIR/contactCounts
inF=$DATADIR/fragmentLists
inB=$DATADIR/biasPerLocus

#####  Settings for EcoRI and HindIII libraries 	  #####
#####  from Duan et al.		   	 		  #####

## meta-fragment size in terms of RE fragments
#howManyFrags=1
## upper and lower bounds on mid-range genomic distances 
distUpThres=250000
distLowThres=10000
## other parameters described in fit-hi-c.py
noOfBins=100
mappabilityThres=1
noOfPasses=1

for i in Duan_yeast_EcoRI Duan_yeast_HindIII; do
   python ../fithic-runner.py -l "$i" -f $inF/$i.gz -i $inI/$i.gz -L $distLowThres -U $distUpThres -b $noOfBins -p $noOfPasses -o outputs/$i --quiet -x interOnly -v
done

#####  Settings for only chromosome 1 of 		#####
#####  human and mouse embryonic stem cell data 	#####
#####  from Dixon et al.	 		  	#####

## meta-fragment size in terms of RE fragments
#howManyFrags=10
## upper and lower bounds on mid-range genomic distances 
distUpThres=5000000
distLowThres=50000
## other parameters described in fit-hi-c.py
noOfBins=200
mappabilityThres=1
noOfPasses=1

# without normalization
for i in Dixon_hESC_HindIII_hg18_combineFrags10_chr1 Dixon_mESC_HindIII_mm9_combineFrags10_chr1; do
   python ../fithic-runner.py -l "$i" -f $inF/$i.gz -i $inI/$i.gz -L $distLowThres -U $distUpThres -b $noOfBins -p $noOfPasses -o outputs/$i --quiet -x intraOnly
done

# without normalization - fixed size windows
noOfBins=50
for i in Dixon_hESC_HindIII_hg18_w40000_chr1; do
   python ../fithic-runner.py -l "$i" -f $inF/$i.gz -i $inI/$i.gz -L $distLowThres -U $distUpThres -b $noOfBins -p $noOfPasses -o outputs/$i --quiet -x interOnly
done

# with normalization - fixed size windows
noOfBins=50
for i in Dixon_hESC_HindIII_hg18_w40000_chr1; do
   fithic -l "$i" -f $inF/$i.gz -i $inI/$i.gz -L $distLowThres -U $distUpThres -b $noOfBins -p $noOfPasses -t $inB/$i.gz -o outputs/$i.afterICE --quiet -x intraOnly
done

echo ""
echo ""
echo "All tests completed successfully. fithic is up and running!"


exit

