#!/bin/bash
# Run CMfinder 0.4 on every job FASTA in a directory (parallel, per-job timeout), then score every
# motif with R-scape (given-structure test; significant covariation at E < 0.05, the default).
# usage: run_cmfinder.sh JOBDIR NPROC TIMEOUT_SEC [JOBLIST]
source ~/mf/bin/activate rs
JD=$1; NP=${2:-16}; TO=${3:-1800}; JL=$4
cd $JD
if [ -n "$JL" ]; then cat $JL; else ls *.fa | sed 's/\.fa$//'; fi | xargs -P $NP -I{} bash -c '
  j={}
  [ -f $j.done ] && exit 0
  mkdir -p $j.cmf && cp $j.fa $j.cmf/ && cd $j.cmf
  timeout '$TO' cmfinder04.pl -skipClustalw -combine -fragmentary -motifList motifs.list -minCandScoreInFinal 10 $j.fa > cmf.log 2>&1
  rc=$?
  for m in $(cat motifs.list 2>/dev/null); do
     mkdir -p rs_$m
     timeout 600 R-scape --outdir rs_$m --nofigures -s $m > rs_$m/out.txt 2>&1 || true
  done
  rm -f *.cand.* *.align.* 2>/dev/null
  echo $rc > ../$j.done
'
echo ALLDONE
