#!/usr/bin/env python
"""Genome-level frame test at CGG sites in marker genes.

For every in-frame CGG in a marker CDS of a sampled genome:
- locate the CDS in the genome (exact match, either strand) and extract the sequence that continues beyond it
- build two chimeric peptides: upstream (up to 40 aa before the CGG, in frame) + downstream 70 aa translated
  (a) in frame 0 from the CGG, or (b) in frame +1 from CGG[1:], i.e. skipping the C
- align both to the marker HMM and record the best domain score; also record the number of stops in the first
  70 downstream codons of each frame
Repeat at random control positions (an in-frame GGN codon) to calibrate. If CGG is a +1 shift site, then:
  score(+1) >> score(0) at CGG, and score(0) > score(+1) at control GGN positions.
usage: gc_frames_genome.py MARKER_ROOT markers.hmm taxmap.tsv genome_paths.tsv NFOCAL NCTRL OUT
"""
import sys, os, glob, gzip, random, numpy as np, pandas as pd, pyhmmer, multiprocessing as mp
ROOT, HMMF, TAX, GP, NF, NC, OUT = sys.argv[1:8]; NF = int(NF); NC = int(NC)
CODE = dict(zip([a + b + c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG'], 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
COMP = str.maketrans('ACGTN', 'TGCAN')
def tr(s): return ''.join(CODE.get(s[i:i + 3], 'X') for i in range(0, len(s) - 2, 3))
def rc(s): return s.translate(COMP)[::-1]
t = pd.read_csv(TAX, sep='\t').set_index('accession')
gp = dict(l.rstrip().split('\t') for l in open(GP))
random.seed(7)
foc = [g for g in t.index[t.family == 'f__CAG-508'] if g[3:] in gp]
ctl = [g for g in t.index[(t['class'] == 'c__Clostridia') & (t.order != 'o__TANB77')] if g[3:] in gp]
foc = random.sample(foc, min(NF, len(foc))); ctl = random.sample(ctl, min(NC, len(ctl)))
markers = sorted(glob.glob(os.path.join(ROOT, 'fna', '*.fna')))

def rfa(fn, keep=None):
    d = {}; n = None; b = []
    op = gzip.open if fn.endswith('.gz') else open
    with op(fn, 'rt') as f:
        for l in f:
            if l[0] == '>':
                if n and (keep is None or n in keep): d[n] = ''.join(b)
                n = l[1:].split()[0]; b = []
            else: b.append(l.strip())
    if n and (keep is None or n in keep): d[n] = ''.join(b)
    return d

_hm = None
def hms():
    global _hm
    if _hm is None:
        _hm = {}
        with pyhmmer.plan7.HMMFile(HMMF) as f:
            for h in f:
                a = h.accession or h.name; a = a.decode() if isinstance(a, bytes) else a; _hm[a.split('.')[0]] = h
    return _hm

MK = None
def work(g):
    hm = hms(); alpha = pyhmmer.easel.Alphabet.amino()
    genome = {k: v.upper() for k, v in rfa(gp[g[3:]]).items()}
    gstr = list(genome.values())
    recs = []; seqs = []
    for mk, cds in MK.get(g, []):
        cds = cds.upper()
        # locate the CDS
        loc = None
        for s in gstr:
            i = s.find(cds)
            if i >= 0: loc = (s, i); break
            j = rc(s).find(cds)
            if j >= 0: loc = (rc(s), j); break
        if loc is None: continue
        s, i0 = loc
        cods = [k for k in range(0, len(cds) - 3, 3)]
        cgg = [k for k in cods if cds[k:k + 3] == 'CGG']
        ggn = [k for k in cods if cds[k:k + 2] == 'GG' and k > 60]
        sites = [(k, 'CGG') for k in cgg] + [(k, 'GGN_ctrl') for k in random.sample(ggn, min(2, len(ggn)))]
        for k, kind in sites:
            gk = i0 + k
            up = tr(s[max(i0, gk - 120):gk])
            d0 = tr(s[gk:gk + 210]); d1 = tr(s[gk + 1:gk + 1 + 210])
            for fr, dn in (('f0', d0), ('f1', d1)):
                pep = (up + dn).replace('*', 'X')
                seqs.append((len(recs), fr, mk, pep))
            recs.append(dict(genome=g, clade='CAG-508' if g in focset else 'other Clostridia', marker=mk, kind=kind,
                             next_nt=s[gk + 3] if gk + 3 < len(s) else '', stops0=d0[:70].count('*'), stops1=d1[:70].count('*'),
                             rel_pos=k / len(cds)))
    # score
    bymk = {}
    for idx, fr, mk, pep in seqs: bymk.setdefault(mk, []).append((idx, fr, pep))
    for mk, lst in bymk.items():
        h = hm.get(mk.split('.')[0]);
        if h is None: continue
        ds = [pyhmmer.easel.TextSequence(name=f'{i}|{fr}'.encode(), sequence=p).digitize(alpha) for i, fr, p in lst]
        blk = pyhmmer.easel.DigitalSequenceBlock(alpha, ds)
        for hits in pyhmmer.hmmsearch([h], blk, E=1e9, domE=1e9, T=-1e9, domT=-1e9, cpus=1):
            for hit in hits:
                nm = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
                i, fr = nm.split('|'); i = int(i)
                recs[i]['score_' + fr] = max(recs[i].get('score_' + fr, -1e9), hit.score)
    return recs

if __name__ == '__main__':
    focset = set(foc); want = set(foc) | set(ctl)
    MK = {}
    for fna in markers:
        mk = os.path.basename(fna)[:-4]
        for g, sq in rfa(fna, want).items(): MK.setdefault(g, []).append((mk, sq))
    with mp.Pool(28) as pool:
        allr = []
        for r in pool.imap_unordered(work, foc + ctl): allr += r
    df = pd.DataFrame(allr); df['delta'] = df.score_f1 - df.score_f0
    df.to_csv(OUT, sep='\t', index=False)
    for (cl, kind), d in df.groupby(['clade', 'kind']):
        print(f'{cl:16s} {kind:9s} n={len(d):6d}  frame+1 better: {(d.delta > 0).mean():.3f}  median delta {d.delta.median():6.1f}  '
              f'stops in 70 codons f0/f1 median {d.stops0.median():.0f}/{d.stops1.median():.0f}  next A/T/G/C '
              + '/'.join(f"{(d.next_nt == x).mean():.2f}" for x in 'ATGC'))
