#!/bin/bash

#DESCRIPTION:
#bash merge-filter_parallelized.sh [inputFile] [resolution] [outputDirectory] [fdr]
#        
#        [inputFile]                Input file of Fit-Hi-C significant interactions
#        [resolution]               Resolution used
#        [outputDirectory]          Output directory to dump output to
#        [fdr]                      False Discovery rate to use when subsetting interactions

inputFile=$1
resolution=$2
outdir=$3
fdr=$4

#JOBHEADER_FILE="~/bin/jobHeader"
UTILITYFOLDER="~/fithic/fithic/utils"
script="$UTILITYFOLDER/CombineNearbyInteraction.py"
mkdir -p $outdir
zcat $inputFile | cut -f1 | sort | uniq > $outdir/chromosomes.used
while read chrom; do
    mkdir -p $outdir/$chrom
    zcat $inputFile | awk '{if(NR!=1){print $0}}'| awk -v c="$chrom" '{if($1==c && $3==c){print $0}}' | awk -v q="$fdr" '{if($7<=q){print $0}}' | gzip > $outdir/$chrom/subset_fithic_"$chrom".gz
    jobFile=$outdir/$chrom/fithic_"$chrom".job 
    #cat $JOBHEADER_FILE > $jobFile
    echo "python3 $script -i $outdir/$chrom/subset_fithic_"$chrom".gz -H 0 -r $resolution -o $outdir/$chrom/postmerged_fithic_"$chrom".gz" >> $jobFile
    #qsub $jobFile
done <$outdir/chromosomes.used
rm $outdir/chromosomes.used
