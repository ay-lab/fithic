interaction=$1
resolution=$2
qvalthresh=$3
outfile=$4

half=$(echo "$resolution/2"|bc)

zcat $interaction | awk -v thresh="${qvalthresh}" -v res="${half}" '{if(NR==1){print "chr1 start1 end1 chr2 start2 end2 cc pval qval"} else { if($7<=thresh) {print $1,$2-half,$2+half,$3,$4-half,$4+half,$5,$6,$7}}}' | gzip > $outfile
