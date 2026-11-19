#!/usr/bin/env python
"""Build CMfinder jobs: gene family (KO) x taxonomic clade.

Each job is a FASTA of 5' leaders (the last UP nt upstream of the start codon plus IN nt into the CDS)
for operon heads that carry the KO, drawn from one GTDB class. Large classes are split by order.
Sequences that are >= 90% identical are collapsed with MMseqs2, and each job is capped at CAP sequences.

usage: make_jobs.py UTR_TSV UTR_FA TAXMAP OUTDIR [KO_LIST]
"""
import sys, os, random, subprocess, collections, tempfile
import pandas as pd

UTSV, UFA, TAX, OUT = sys.argv[1:5]
KOS = set(l.strip() for l in open(sys.argv[5])) if len(sys.argv) > 5 else None
UP, IN, MINUP, MINN, CAP, BIG = 400, 20, 80, 8, 60, 120
random.seed(1)

def readfa(fn):
    d = {}; n = None
    for l in open(fn):
        if l[0] == '>': n = l[1:].strip(); d[n] = []
        else: d[n].append(l.strip())
    return {k: ''.join(v) for k, v in d.items()}

u = pd.read_csv(UTSV, sep='\t')
u = u[u.up_len >= MINUP]
tax = pd.read_csv(TAX, sep='\t').set_index('accession')
u['acc'] = u.genome
tx = tax.reindex(u.acc)
u['cls'] = tx['class'].values; u['ord'] = tx['order'].values; u['phy'] = tx['phylum'].values
seqs = readfa(UFA)

def trim(uid, up_len):
    s = seqs[uid]; cds = len(s) - up_len
    s = s[max(0, up_len - UP): up_len + min(IN, cds)]
    return s

def derep(recs):
    if len(recs) < 3: return recs
    with tempfile.TemporaryDirectory() as td:
        fa = os.path.join(td, 'in.fa')
        with open(fa, 'w') as f:
            for k, s in recs: f.write(f'>{k}\n{s}\n')
        subprocess.run(['mmseqs', 'easy-cluster', fa, os.path.join(td, 'c'), os.path.join(td, 't'),
                        '--min-seq-id', '0.9', '-c', '0.8', '--cov-mode', '1', '--dbtype', '2', '-v', '0',
                        '--threads', '1'], check=True)
        reps = set(l.split('\t')[0] for l in open(os.path.join(td, 'c_cluster.tsv')))
    return [r for r in recs if r[0] in reps]

rows = []
u['kolist'] = u.kos.str.split(',')
ex = u.explode('kolist')
if KOS: ex = ex[ex.kolist.isin(KOS)]
os.makedirs(OUT, exist_ok=True)
for ko, g in ex.groupby('kolist'):
    for cls, gc in g.groupby('cls'):
        parts = [(cls, gc)] if len(gc) <= BIG else list(gc.groupby('ord'))
        leftovers = []
        for name, gp in parts:
            if len(gp) < MINN:
                leftovers.append(gp); continue
            rows.append((ko, name, gp))
        if leftovers:
            lo = pd.concat(leftovers)
            if len(lo) >= MINN: rows.append((ko, cls + '_other', lo))
n = 0
manifest = open(os.path.join(OUT, 'manifest.tsv'), 'w')
manifest.write('job\tko\tclade\tn_raw\tn_derep\tn_used\n')
for ko, clade, gp in rows:
    recs = [(r.uid, trim(r.uid, r.up_len)) for r in gp.itertuples()]
    dr = derep(recs)
    if len(dr) < MINN: continue
    use = random.sample(dr, CAP) if len(dr) > CAP else dr
    job = f'{ko}__{clade.replace(";", "_")}'
    with open(os.path.join(OUT, job + '.fa'), 'w') as f:
        for k, s in use: f.write(f'>{k}\n{s}\n')
    manifest.write(f'{job}\t{ko}\t{clade}\t{len(recs)}\t{len(dr)}\t{len(use)}\n'); n += 1
print('jobs', n)
