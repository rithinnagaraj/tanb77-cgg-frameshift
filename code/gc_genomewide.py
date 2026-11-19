#!/usr/bin/env python
"""Genome-wide extent of CGG-associated +1 frameshifting.
For each genome (Prodigal CDS from GTDB faa headers plus the genome sequence):
- n_genes
- n_cgg: number of in-frame CGG codons
- n_cgg_term: CDSs whose last 15 codons (before the stop) contain a CGG
- n_shift_orf: of those, how many have a +1-frame ORF from CGG[1:] that runs >= 60 codons beyond the annotated
  stop (i.e., translation would continue in +1)
- control: the same metrics for AGG (an Arg codon that does have a tRNA) and for GGG
usage: gc_genomewide.py LIST.tsv(genome<TAB>fna.gz<TAB>faa.gz) OUT NPROC
"""
import sys, gzip, multiprocessing as mp
CODE = dict(zip([a + b + c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG'], 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
COMP = str.maketrans('ACGTN', 'TGCAN')
def rfa(fn):
    d = {}; n = None; b = []
    with gzip.open(fn, 'rt') as f:
        for l in f:
            if l[0] == '>':
                if n: d[n] = ''.join(b)
                n = l[1:].rstrip(); b = []
            else: b.append(l.strip())
    if n: d[n] = ''.join(b)
    return d
def orflen(s):
    for i in range(0, len(s) - 2, 3):
        if CODE.get(s[i:i + 3]) == '*': return i // 3
    return len(s) // 3
def work(row):
    gid, fna, faa = row
    g = {k.split()[0]: v.upper() for k, v in rfa(fna).items()}
    res = {'genes': 0}
    for cod in ('CGG', 'AGG', 'GGG'):
        for k in ('n', 'term', 'shift'): res[f'{cod}_{k}'] = 0
    for h in rfa(faa):
        p = h.split(' # '); pid = p[0].strip(); s, e, st = int(p[1]), int(p[2]), int(p[3])
        c = pid.rsplit('_', 1)[0]
        if c not in g: continue
        seq = g[c]
        if st == 1:
            cds = seq[s - 1:e]; ext = seq[s - 1:min(len(seq), e + 600)]
        else:
            cds = seq[s - 1:e].translate(COMP)[::-1]; ext = seq[max(0, s - 1 - 600):e].translate(COMP)[::-1]
        if 'partial=00' not in h and 'partial' in h: pass
        res['genes'] += 1
        L = len(cds) - 3
        for cod in ('CGG', 'AGG', 'GGG'):
            pos = [k for k in range(0, L, 3) if cds[k:k + 3] == cod]
            res[f'{cod}_n'] += len(pos)
            tp = [k for k in pos if k >= L - 45]
            if tp:
                res[f'{cod}_term'] += 1
                k = tp[-1]
                if orflen(ext[k + 1:]) * 3 + k + 1 >= len(cds) + 180: res[f'{cod}_shift'] += 1
    return gid + '\t' + '\t'.join(str(res[k]) for k in sorted(res))
if __name__ == '__main__':
    rows = [l.rstrip('\n').split('\t') for l in open(sys.argv[1])]
    keys = sorted(['genes'] + [f'{c}_{k}' for c in ('CGG', 'AGG', 'GGG') for k in ('n', 'term', 'shift')])
    with mp.Pool(int(sys.argv[3])) as pool, open(sys.argv[2], 'w') as fo:
        fo.write('genome\t' + '\t'.join(keys) + '\n')
        for r in pool.imap_unordered(work, rows): fo.write(r + '\n')
