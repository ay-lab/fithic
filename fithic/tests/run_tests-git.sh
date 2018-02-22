#!/bin/bash -ex

set -o nounset
set -o pipefail
set -o errexit 

BINDIR=../bin
DATADIR=data
inI=$DATADIR/contactCounts
inF=$DATADIR/fragmentLists
inB=$DATADIR/biasPerLocus


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
   python ../fithic-runner.py -l "$i" -f $inF/$i.gz -i $inI/$i.gz -o outputs/$i -L $distLowThres -U $distUpThres -b $noOfBins -p $noOfPasses --quiet -x intraOnly
done

# without normalization - fixed size windows
noOfBins=50
for i in Dixon_hESC_HindIII_hg18_w40000_chr1; do
   python ../fithic-runner.py -l "$i" -f $inF/$i.gz -i $inI/$i.gz -L $distLowThres -U $distUpThres -b $noOfBins -p $noOfPasses -o outputs/$i --quiet -x All
done

# with normalization - fixed size windows
noOfBins=50
for i in Dixon_hESC_HindIII_hg18_w40000_chr1; do
   python ../fithic-runner.py -l "$i" -f $inF/$i.gz -i $inI/$i.gz -L $distLowThres -U $distUpThres -b $noOfBins -p $noOfPasses -t $inB/$i.gz -o outputs/$i.afterICE --quiet -x intraOnly
done

noOfBins=200
#without normalization - interOnly and All
for i in hIMR90_100kb_hg19; do
    python ../fithic-runner.py -l "$i" -i $inI/${i}_cc.gz -f $inF/${i}_frags.gz -b $noOfBins -p $noOfPasses -o outputs/${i}.interOnly -x interOnly
    python ../fithic-runner.py -l "$i" -i $inI/${i}_cc.gz -f $inF/${i}_frags.gz -b $noOfBins -p $noOfPasses -o outputs/${i}.all -x All
done

#USE BELOW IF WANTING TO TRY NONHUMAN/MOUSE ORGANISMS
#for i in RINGS; do
#    python ../fithic-runner.py -l "$i" -i $inI/${i}.gz -f $inF/${i}.gz -b $noOfBins -p $noOfPasses -o outputs/${i} -x All
#done

echo ""
echo ""
echo "All tests completed successfully. fithic is up and running!"


exit

