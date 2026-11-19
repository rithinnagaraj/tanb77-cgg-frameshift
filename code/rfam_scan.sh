#!/bin/bash
# Annotate every extracted 5' leader with Rfam 15 (gathering thresholds) and also keep sub-threshold
# hits (E < 1e-3) so that "variants of known classes" at blind-spot genes can be inspected.
# usage: rfam_scan.sh UTR_FASTA OUTPREFIX NCHUNK
source ~/mf/bin/activate rs
FA=$1; OUT=$2; N=${3:-16}
RF=~/data/rfam
mkdir -p ${OUT}_chunks
seqkit split2 -p $N -O ${OUT}_chunks -f $FA > /dev/null 2>&1
ls ${OUT}_chunks/*.fa* | xargs -P $N -I{} sh -c 'cmscan --cpu 1 --rfam --nohmmonly --fmt 2 --clanin '$RF'/Rfam.clanin --oclan -E 1e-3 --tblout {}.tbl '$RF'/Rfam.cm {} > /dev/null'
cat ${OUT}_chunks/*.tbl | grep -v "^#" > ${OUT}.tbl
echo done
