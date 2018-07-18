#!/bin/bash -ex

set -o nounset
set -o pipefail
set -o errexit 

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
    fithic -r 0 -l "$i" -f $inF/$i.gz -i $inI/$i.gz -o outputs/$i -L $distLowThres -U $distUpThres -b $noOfBins -p $noOfPasses -x intraOnly
done

# without normalization - fixed size windows
noOfBins=50
for i in Dixon_hESC_HindIII_hg18_w40000_chr1; do
    fithic -r 40000 -l "$i" -f $inF/$i.gz -i $inI/$i.gz -L $distLowThres -U $distUpThres -b $noOfBins -p $noOfPasses -o outputs/$i -x All
done

# with normalization - fixed size windows
noOfBins=50
for i in Dixon_hESC_HindIII_hg18_w40000_chr1; do
    fithic -r 40000 -l "$i" -f $inF/$i.gz -i $inI/$i.gz -L $distLowThres -U $distUpThres -b $noOfBins -p $noOfPasses -t $inB/$i.gz -o outputs/$i.afterICE -x intraOnly
done

noOfBins=200
# with normalization - interOnly and All
for i in Dixon_IMR90_HindIII_hg19_w100000; do
    fithic -r 100000 -l "$i" -i $inI/$i.gz -f $inF/$i.gz -t $inB/$i.gz -b $noOfBins -p $noOfPasses -o outputs/${i}.interOnly -x interOnly
    fithic -r 100000 -l "$i" -i $inI/$i.gz -f $inF/$i.gz -t $inB/$i.gz -b $noOfBins -p $noOfPasses -o outputs/${i}.all -x All
done

# without normalization - All, nonhuman
for i in Ay_Rings_MboI_Pfal_w10000; do
    fithic -r 10000 -l "$i" -i $inI/$i.gz -f $inF/$i.gz -b $noOfBins -p $noOfPasses -o outputs/${i} -x All
done

echo ""
echo ""
echo "All tests completed successfully. Fit-Hi-C is up and running!"


exit

