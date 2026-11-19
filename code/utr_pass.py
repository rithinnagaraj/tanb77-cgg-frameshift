#!/usr/bin/env python
"""Pass 2: cut the 5' leader (up to L nt upstream of the operon head's start codon, bounded by the
upstream neighbour gene, plus D nt into the CDS) for every operon head that pass 1 found.

usage: utr_pass.py HITS_TSV GENOME_PATHS_TSV OUT_FASTA OUT_TSV [L] [D] [NPROC]
GENOME_PATHS_TSV: accession<TAB>path_to_fna.gz
"""
import sys, gzip, multiprocessing as mp
import pandas as pd

HITS, GP, OFA, OTSV = sys.argv[1:5]
L = int(sys.argv[5]) if len(sys.argv) > 5 else 600
D = int(sys.argv[6]) if len(sys.argv) > 6 else 60
NP = int(sys.argv[7]) if len(sys.argv) > 7 else 16
COMP = str.maketrans('ACGTNacgtnRYKMSWBDHV', 'TGCANtgcanYRMKSWVHDB')

def rc(s): return s.translate(COMP)[::-1]

def load_fasta(fn):
    seqs = {}; name = None; buf = []
    with gzip.open(fn, 'rt') as fh:
        for line in fh:
            if line[0] == '>':
                if name: seqs[name] = ''.join(buf)
                name = line[1:].split()[0]; buf = []
            else: buf.append(line.strip())
    if name: seqs[name] = ''.join(buf)
    return seqs

def work(args):
    gid, path, rows = args
    try: seqs = load_fasta(path)
    except Exception as ex: return [], [f'#ERR {gid} {ex}']
    fa, tab = [], []
    for (contig, atg, st, ub, head, kos) in rows:
        s = seqs.get(contig)
        if s is None: continue
        n = len(s)
        if st == 1:
            lo = max(1, atg - L, (ub + 1) if ub > 0 else 1)
            hi = min(n, atg + D - 1)
            sub = s[lo - 1:hi]; igr = atg - ub - 1 if ub > 0 else atg - 1
            edge = ub <= 0
        else:
            hi = min(n, atg + L, (ub - 1) if ub > 0 else n)
            lo = max(1, atg - D + 1)
            sub = rc(s[lo - 1:hi]); igr = ub - atg - 1 if ub > 0 else n - atg
            edge = ub <= 0
        up_len = len(sub) - (min(D, atg) if st == 1 else min(D, n - atg + 1))
        if up_len < 20: continue
        uid = f'{gid}|{head}'
        fa.append(f'>{uid}\n{sub.upper()}\n')
        tab.append(f'{uid}\t{gid}\t{contig}\t{st}\t{atg}\t{igr}\t{up_len}\t{int(edge)}\t{kos}')
    return fa, tab

if __name__ == '__main__':
    h = pd.read_csv(HITS, sep='\t', comment='#')
    gp = dict(l.rstrip('\n').split('\t') for l in open(GP))
    g = h.groupby(['genome', 'head']).agg(contig=('contig', 'first'), atg=('head_atg', 'first'),
            st=('head_strand', 'first'), ub=('up_bound', 'first'),
            kos=('ko', lambda x: ','.join(sorted(set(x))))).reset_index()
    jobs = []
    for gid, sub in g.groupby('genome'):
        acc = gid[3:] if gid[:3] in ('RS_', 'GB_') else gid
        if acc not in gp: continue
        jobs.append((gid, gp[acc], [(r.contig, int(r.atg), int(r.st), int(r.ub), r.head, r.kos) for r in sub.itertuples()]))
    print('genomes', len(jobs), 'heads', len(g), flush=True)
    with mp.Pool(NP) as pool, open(OFA, 'w') as fa, open(OTSV, 'w') as tb:
        tb.write('uid\tgenome\tcontig\tstrand\tatg\tigr_len\tup_len\tcontig_edge\tkos\n')
        for i, (f, t) in enumerate(pool.imap_unordered(work, jobs, chunksize=8)):
            fa.writelines(f); tb.write(''.join(x + '\n' for x in t if not x.startswith('#')))
            if i % 2000 == 0: print(i, flush=True)
