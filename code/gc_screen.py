#!/usr/bin/env python
"""Fast, Codetta-like genetic-code screen over the GTDB marker genes (bac120 + ar53).

For every species-representative genome:
- Each marker protein (as GTDB translated it) is aligned to its profile HMM (Pfam or TIGRFAM) with pyhmmer.
- Every residue in a match state k is mapped back to its codon in the marker CDS.
- For every (genome, codon) the script accumulates
      S[a]     = sum_k log e_k(a)   over the aligned match columns (a = 20 amino acids)
      N        = number of aligned match columns
      C[a]     = number of those columns whose consensus (max emission >= 0.5) is amino acid a
Workers split the work by marker; each keeps full-size accumulators, and the main process sums them.

usage: gc_screen.py MARKER_ROOTS(comma) HMMDB(comma) OUT.npz NPROC
"""
import sys, os, glob, numpy as np, multiprocessing as mp
import pyhmmer

ROOTS = sys.argv[1].split(','); HMMDBS = sys.argv[2].split(','); OUT = sys.argv[3]; NP = int(sys.argv[4])
BASES = 'TCAG'
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
CIDX = {c: i for i, c in enumerate(CODONS)}

def read_fa(fn):
    seqs = {}; n = None; buf = []
    for l in open(fn):
        if l[0] == '>':
            if n: seqs[n] = ''.join(buf)
            n = l[1:].split()[0]; buf = []
        else: buf.append(l.strip())
    if n: seqs[n] = ''.join(buf)
    return seqs

def marker_list():
    ms = []
    for r in ROOTS:
        for f in sorted(glob.glob(os.path.join(r, 'faa', '*.faa'))):
            m = os.path.basename(f)[:-4]
            ms.append((m, f, os.path.join(r, 'fna', m + '.fna')))
    return ms

def load_hmm(acc):
    for db in HMMDBS:
        with pyhmmer.plan7.HMMFile(db) as f:
            for h in f:
                a = h.accession or h.name
                a = a.decode() if isinstance(a, bytes) else a
                if a.split('.')[0] == acc: return h
    return None

GIDX = None
def init(gidx_file):
    global GIDX
    GIDX = {g: i for i, g in enumerate(open(gidx_file).read().split())}

def worker(markers):
    G = len(GIDX)
    S = np.zeros((G, 64, 20), dtype=np.float32)
    N = np.zeros((G, 64), dtype=np.int32)
    C = np.zeros((G, 64, 20), dtype=np.int16)
    alpha = pyhmmer.easel.Alphabet.amino()
    aligner = pyhmmer.plan7.TraceAligner()
    log = []
    for m, faa, fna in markers:
        h = load_hmm(m.split('.')[0])
        if h is None: log.append(f'{m} NOHMM'); continue
        em = np.clip(np.array(h.match_emissions)[1:, :20], 1e-6, 1)
        lem = np.log(em).astype(np.float32); cons = em.argmax(1); consok = em.max(1) >= 0.5
        P = read_fa(faa); Nn = read_fa(fna)
        names = [g for g in P if g in Nn and g in GIDX]
        for b in range(0, len(names), 5000):
            sub = names[b:b + 5000]
            seqs = [pyhmmer.easel.TextSequence(name=g.encode(), sequence=P[g].replace('*', '').replace('X', 'A')).digitize(alpha) for g in sub]
            blk = pyhmmer.easel.DigitalSequenceBlock(alpha, seqs)
            traces = aligner.compute_traces(h, blk)
            msa = aligner.align_traces(h, blk, traces, all_consensus_cols=True, trim=False)
            ref = msa.reference
            isx = [ch == 'x' for ch in ref]
            midx = []; k = 0
            for x in isx:
                midx.append(k if x else -1)
                if x: k += 1
            for g, row in zip(sub, msa.alignment):
                nt = Nn[g].upper(); gi = GIDX[g]
                ri = -1
                for col, ch in enumerate(row):
                    if not ch.isalpha(): continue
                    ri += 1
                    if not (isx[col] and ch.isupper()): continue
                    kk = midx[col]
                    c = CIDX.get(nt[3 * ri: 3 * ri + 3])
                    if c is None: continue
                    if consok[kk]:
                        S[gi, c] += lem[kk]; N[gi, c] += 1; C[gi, c, cons[kk]] += 1
        log.append(f'{m} n={len(names)}')
        print(log[-1], flush=True)
    return S, N, C

if __name__ == '__main__':
    ms = marker_list()
    gidx_file = OUT + '.genomes.txt'
    if not os.path.exists(gidx_file):
        gs = set()
        for m, faa, fna in ms:
            gs |= set(l[1:].split()[0] for l in open(faa) if l[0] == '>')
        open(gidx_file, 'w').write('\n'.join(sorted(gs)) + '\n')
    groups = [ms[i::NP] for i in range(NP)]
    with mp.Pool(NP, initializer=init, initargs=(gidx_file,)) as pool:
        tot = None
        for S, N, C in pool.imap_unordered(worker, groups):
            if tot is None: tot = [S, N, C.astype(np.int32)]
            else: tot[0] += S; tot[1] += N; tot[2] += C
    np.savez_compressed(OUT, S=tot[0], N=tot[1], C=tot[2])
    print('saved', OUT)
