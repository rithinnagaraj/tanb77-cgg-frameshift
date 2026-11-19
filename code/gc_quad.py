#!/usr/bin/env python
"""Test the CGGN quadruplet / +1 frameshift hypothesis in a focal clade versus a control clade.

For every in-frame CGG codon in the marker CDSs (GTDB nucleotide marker genes) of focal and control genomes:
- the next nucleotide N (in CGGN)
- the distance (in codons) to the next stop codon in frame 0 (the annotated frame) and in frame +1
  (starting right after the extra C, i.e. reading GGN...)
- homology test: translate 15 codons downstream in frame 0 and in frame +1, and score both against the same
  marker's reference consensus. A pyhmmer local alignment of the downstream peptide to the marker HMM is the
  proxy: take the best domain score of each peptide.
Reports summary tables per clade.

usage: gc_quad.py MARKER_ROOT markers.hmm taxmap.tsv OUTPREFIX
"""
import sys, glob, os, collections, numpy as np, pandas as pd, pyhmmer
ROOT, HMMF, TAX, OUT = sys.argv[1:5]
t = pd.read_csv(TAX, sep='\t').set_index('accession')
focal = set(t.index[t.family == 'f__CAG-508'])
ctrl_pool = list(t.index[(t['class'] == 'c__Clostridia') & (t.order != 'o__TANB77')])
rng = np.random.default_rng(0); ctrl = set(rng.choice(ctrl_pool, 2000, replace=False))
CODE = {a + b + c: aa for (a, b, c), aa in zip([(x, y, z) for x in 'TCAG' for y in 'TCAG' for z in 'TCAG'],
        'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')}
def tr(s): return ''.join(CODE.get(s[i:i + 3], 'X') for i in range(0, len(s) - 2, 3))
def stopdist(s):
    p = tr(s); i = p.find('*'); return i if i >= 0 else len(p)
def rfa(fn):
    d = {}; n = None; b = []
    for l in open(fn):
        if l[0] == '>':
            if n: d[n] = ''.join(b)
            n = l[1:].split()[0]; b = []
        else: b.append(l.strip())
    if n: d[n] = ''.join(b)
    return d
hm = {}
with pyhmmer.plan7.HMMFile(HMMF) as f:
    for h in f:
        a = h.accession or h.name; a = a.decode() if isinstance(a, bytes) else a; hm[a.split('.')[0]] = h
alpha = pyhmmer.easel.Alphabet.amino()
rows = []
for fna in sorted(glob.glob(os.path.join(ROOT, 'fna', '*.fna'))):
    mk = os.path.basename(fna)[:-4]; h = hm.get(mk.split('.')[0])
    N = rfa(fna)
    peps = []
    for g, nt in N.items():
        if g not in focal and g not in ctrl: continue
        nt = nt.upper()
        for k in range(0, len(nt) - 3, 3):
            if nt[k:k + 3] != 'CGG': continue
            down0 = nt[k + 3:k + 3 + 120]; down1 = nt[k + 1:k + 1 + 120]
            rec = dict(genome=g, clade='CAG-508' if g in focal else 'other Clostridia', marker=mk, codon_idx=k // 3,
                       rel_pos=round(k / len(nt), 3), next_nt=nt[k + 3] if k + 3 < len(nt) else '',
                       stop0=stopdist(down0), stop1=stopdist(down1), pep0=tr(nt[k + 3:k + 3 + 45]), pep1=tr(nt[k + 1:k + 1 + 45]))
            rows.append(rec); peps.append(len(rows) - 1)
    if h is not None and peps:
        # score the frame-0 vs frame+1 downstream peptides against the marker HMM
        seqs = []
        for i in peps:
            for fr in ('pep0', 'pep1'):
                p = rows[i][fr].replace('*', 'X')
                if len(p) < 5: p = p + 'A' * (5 - len(p))
                seqs.append(pyhmmer.easel.TextSequence(name=f'{i}|{fr}'.encode(), sequence=p).digitize(alpha))
        blk = pyhmmer.easel.DigitalSequenceBlock(alpha, seqs)
        best = {}
        for hits in pyhmmer.hmmsearch([h], blk, E=1e9, domE=1e9, T=-1e9, domT=-1e9):
            for hit in hits:
                nm = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
                best[nm] = max(best.get(nm, -1e9), hit.score)
        for i in peps:
            rows[i]['score0'] = best.get(f'{i}|pep0', np.nan); rows[i]['score1'] = best.get(f'{i}|pep1', np.nan)
    print(mk, len(rows), flush=True)
df = pd.DataFrame(rows)
df.to_csv(OUT + '_cgg_frames.tsv', sep='\t', index=False)
for cl, d in df.groupby('clade'):
    print(f'\n== {cl}: CGG occurrences {len(d)} in {d.genome.nunique()} genomes')
    print(' next nt:', d.next_nt.value_counts(normalize=True).round(3).to_dict())
    print(' median codons to stop: frame0 %.0f  frame+1 %.0f' % (d.stop0.median(), d.stop1.median()))
    print(' frame+1 scores higher than frame0 (HMM):', round(float((d.score1 > d.score0).mean()), 3),
          ' median score0 %.1f score1 %.1f' % (d.score0.median(), d.score1.median()))
    print(' within last 5% of CDS (annotated end):', round(float((d.rel_pos > 0.95).mean()), 3))
