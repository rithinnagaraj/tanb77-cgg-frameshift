#!/usr/bin/env python
"""Site-level evidence for one codon in a focal clade.

Every occurrence of CODON in the marker genes of focal genomes is aligned to the marker HMM, and its column
is recorded. For the same marker columns, the amino acid used by reference (non-focal) genomes is tallied:
how conserved is Gly or Arg at the sites where the focal clade puts CODON?

usage: gc_sites.py MARKER_ROOT HMM.hmm taxmap.tsv CODON FOCAL_RANK FOCAL_NAME REF_RANK REF_NAME OUTPREFIX
"""
import sys, os, glob, numpy as np, pandas as pd, pyhmmer, collections
ROOT, HMMF, TAX, CODON, FR, FN, RR, RN, OUT = sys.argv[1:10]
tax = pd.read_csv(TAX, sep='\t').set_index('accession')
focal = set(tax.index[tax[FR] == FN]); ref = set(tax.index[(tax[RR] == RN) & (tax[FR] != FN)])
hmms = {}
with pyhmmer.plan7.HMMFile(HMMF) as f:
    for h in f:
        a = h.accession or h.name; a = a.decode() if isinstance(a, bytes) else a
        hmms[a.split('.')[0]] = h
alpha = pyhmmer.easel.Alphabet.amino(); al = pyhmmer.plan7.TraceAligner()
AA = 'ACDEFGHIKLMNPQRSTVWY'

def rfa(fn):
    d = {}; n = None; b = []
    for l in open(fn):
        if l[0] == '>':
            if n: d[n] = ''.join(b)
            n = l[1:].split()[0]; b = []
        else: b.append(l.strip())
    if n: d[n] = ''.join(b)
    return d

occ = []; refcols = collections.defaultdict(collections.Counter)
for faa in sorted(glob.glob(os.path.join(ROOT, 'faa', '*.faa'))):
    m = os.path.basename(faa)[:-4]; h = hmms.get(m.split('.')[0])
    if h is None: continue
    em = np.array(h.match_emissions)[1:, :20]
    P = rfa(faa); N = rfa(os.path.join(ROOT, 'fna', m + '.fna'))
    names = [g for g in P if g in N and (g in focal or g in ref)]
    # cap reference genomes for speed
    fnames = [g for g in names if g in focal]; rnames = [g for g in names if g in ref][:3000]
    use = fnames + rnames
    if not fnames: continue
    seqs = [pyhmmer.easel.TextSequence(name=g.encode(), sequence=P[g].replace('*', '').replace('X', 'A')).digitize(alpha) for g in use]
    blk = pyhmmer.easel.DigitalSequenceBlock(alpha, seqs)
    msa = al.align_traces(h, blk, al.compute_traces(h, blk), all_consensus_cols=True, trim=False)
    isx = [c == 'x' for c in msa.reference]; midx = []; k = 0
    for x in isx: midx.append(k if x else -1); k += x
    for g, row in zip(use, msa.alignment):
        nt = N[g].upper(); ri = -1
        for col, ch in enumerate(row):
            if not ch.isalpha(): continue
            ri += 1
            if not (isx[col] and ch.isupper()): continue
            kk = midx[col]
            if g in ref: refcols[(m, kk)][ch] += 1
            elif nt[3 * ri:3 * ri + 3] == CODON:
                occ.append(dict(genome=g, genus=tax.at[g, 'genus'], marker=m, col=kk, aa_translated=ch,
                                hmm_consensus=AA[em[kk].argmax()], hmm_pG=round(float(em[kk][AA.index('G')]), 3),
                                hmm_pR=round(float(em[kk][AA.index('R')]), 3)))
    print(m, len(fnames), len(occ), flush=True)
oc = pd.DataFrame(occ)
def refstats(r):
    c = refcols.get((r.marker, r.col), collections.Counter()); n = sum(c.values())
    return pd.Series(dict(ref_n=n, ref_fracG=c['G'] / n if n else np.nan, ref_fracR=c['R'] / n if n else np.nan,
                          ref_top=c.most_common(1)[0][0] if n else ''))
oc = pd.concat([oc, oc.apply(refstats, axis=1)], axis=1)
oc.to_csv(OUT + '_sites.tsv', sep='\t', index=False)
print('occurrences', len(oc))
print('reference residue at focal CGG sites: Gly-majority', (oc.ref_top == 'G').sum(), 'Arg-majority', (oc.ref_top == 'R').sum(), 'other', ((oc.ref_top != 'G') & (oc.ref_top != 'R')).sum())
print('sites where reference Gly >= 95%:', (oc.ref_fracG >= 0.95).sum(), ' Arg >= 95%:', (oc.ref_fracR >= 0.95).sum())
