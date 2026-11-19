#!/usr/bin/env python
"""Pass 1: search every GTDB representative proteome with the selected KOfam HMMs.

For each hit gene the script records the operon head (walking upstream on the same strand while the
intergenic gap is < OPGAP nt), the head's start coordinate and the boundary of the upstream neighbour.
That gives pass 2 everything it needs to cut the 5' leader from the genome.

usage: hmm_pass.py HMMFILE THRESH_TSV FAA_DIR OUT_TSV [NPROC] [LIMIT]
"""
import sys, os, gzip, glob, multiprocessing as mp
import pyhmmer

OPGAP = 50          # max same-strand gap to call two genes co-operonic
HMMFILE, THR, FAADIR, OUT = sys.argv[1:5]
NPROC = int(sys.argv[5]) if len(sys.argv) > 5 else 16
LIMIT = int(sys.argv[6]) if len(sys.argv) > 6 else 0

_hmms = None; _thr = None; _alpha = pyhmmer.easel.Alphabet.amino()

def init():
    global _hmms, _thr
    with pyhmmer.plan7.HMMFile(HMMFILE) as f:
        _hmms = list(f)
    _thr = {}
    for line in open(THR):
        ko, t, st = line.rstrip('\n').split('\t')[:3]
        _thr[ko] = (float(t), st)

def parse(fn):
    genes = []   # (contig, start, end, strand, protid, seq)
    name = None; seq = []
    with gzip.open(fn, 'rt') as fh:
        for line in fh:
            if line[0] == '>':
                if name: genes.append(name + (''.join(seq),))
                p = line[1:].split(' # ')
                pid = p[0].strip(); s, e, st = int(p[1]), int(p[2]), int(p[3])
                contig = pid.rsplit('_', 1)[0]
                name = (contig, s, e, st, pid); seq = []
            else:
                seq.append(line.strip().rstrip('*'))
        if name: genes.append(name + (''.join(seq),))
    return genes

def work(fn):
    gid = os.path.basename(fn).replace('_protein.faa.gz', '')
    try:
        genes = parse(fn)
    except Exception as ex:
        return [f'#ERR\t{gid}\t{ex}']
    idx = {g[4]: i for i, g in enumerate(genes)}
    seqs = [pyhmmer.easel.TextSequence(name=g[4].encode(), sequence=g[5].replace('*', 'X')).digitize(_alpha)
            for g in genes]
    blk = pyhmmer.easel.DigitalSequenceBlock(_alpha, seqs)
    best = {}  # protid -> (score, ko)
    for hits in pyhmmer.hmmsearch(_hmms, blk, cpus=1, E=1e-5):
        ko = hits.query.name.decode() if isinstance(hits.query.name, bytes) else hits.query.name
        t, st = _thr.get(ko, (0.0, 'full'))
        for h in hits:
            sc = h.score if st == 'full' else max(d.score for d in h.domains)
            if sc < t: continue
            pid = h.name.decode() if isinstance(h.name, bytes) else h.name
            if pid not in best or sc > best[pid][0]:
                best[pid] = (sc, ko)
    out = []
    n = len(genes)
    for pid, (sc, ko) in best.items():
        i = idx[pid]; c, s, e, st = genes[i][:4]
        # walk to operon head
        j = i; nop = 0
        while True:
            k = j - st                     # upstream neighbour index
            if k < 0 or k >= n or genes[k][0] != c or genes[k][3] != st: break
            gap = (genes[j][1] - genes[k][2]) if st == 1 else (genes[k][1] - genes[j][2])
            if gap >= OPGAP: break
            j = k; nop += 1
        hc, hs, he, hst, hp = genes[j][:5]
        # upstream boundary: nearest gene on either strand upstream of head (genes sorted by start per contig)
        k = j - hst
        if 0 <= k < n and genes[k][0] == hc:
            ub = genes[k][2] if hst == 1 else genes[k][1]
            ust = genes[k][3]
        else:
            ub = -1; ust = 0
        tss = hs if hst == 1 else he          # start codon coordinate (1-based)
        out.append(f'{gid}\t{ko}\t{pid}\t{sc:.1f}\t{c}\t{s}\t{e}\t{st}\t{hp}\t{tss}\t{hst}\t{ub}\t{ust}\t{nop}')
    return out

if __name__ == '__main__':
    files = sorted(glob.glob(os.path.join(FAADIR, '**', '*_protein.faa.gz'), recursive=True))
    if LIMIT: files = files[:LIMIT]
    with mp.Pool(NPROC, initializer=init) as pool, open(OUT, 'w') as fo:
        fo.write('genome\tko\tprot\tscore\tcontig\tstart\tend\tstrand\thead\thead_atg\thead_strand\tup_bound\tup_strand\tn_up_in_operon\n')
        for i, res in enumerate(pool.imap_unordered(work, files, chunksize=4)):
            for r in res: fo.write(r + '\n')
            if i % 5000 == 0:
                print(i, len(files), flush=True)
