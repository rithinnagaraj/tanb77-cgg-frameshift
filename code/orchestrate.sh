#!/bin/bash
# After pass 1: extract the leaders, build the jobs, run CMfinder + R-scape (targets first, then
# positive controls, then nulls), and parse the results.
source ~/mf/bin/activate rs
cd ~/work
until grep -q DONE pass1.log; do sleep 60; done
mkdir -p full
python utr_pass.py pass1_genus.tsv genome_paths_all.tsv full/utr.fa full/utr.tsv 600 60 16 > full/utr.log 2>&1
python make_jobs.py full/utr.tsv full/utr.fa taxmap.tsv full/jobs > full/jobs.log 2>&1
python - <<'PY'
import pandas as pd
f = pd.read_csv('families.tsv', sep='\t')
grp = dict(zip(f.ko, f.group))
m = pd.read_csv('full/jobs/manifest.tsv', sep='\t')
def rank(ko):
    g = grp.get(ko, 'NULL')
    return 0 if not g.startswith('POS') and g != 'NULL' else (1 if g.startswith('POS') else 2)
m['r'] = m.ko.map(rank)
m = m.sort_values(['r', 'n_used'], ascending=[True, False])
m.job.to_csv('full/joblist.txt', index=False, header=False)
print(m.r.value_counts())
PY
bash run_cmfinder.sh full/jobs 16 1500 ~/work/full/joblist.txt > full/cmf.log 2>&1
python parse_motifs.py full/jobs full/utr.tsv full/motifs.tsv full/motifs.fa > full/parse.log 2>&1
echo ORCH_DONE >> full/cmf.log
