HIC_DUMPED=$1
CHR1=$2
CHR2=$3
OUT=$4

cat $HIC_DUMPED | awk -v chr1=$CHR1 -v chr2=$CHR2 '{ printf "%s\t%s\t%s\t%s\t%s\n", chr1,$1,chr2,$2,$3}' | gzip > $OUT
