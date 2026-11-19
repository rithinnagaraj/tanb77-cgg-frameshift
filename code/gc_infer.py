#!/usr/bin/env python
"""Infer codon meanings per genome and per clade from the marker screen.

Per genome and codon:
- posterior P(a) is proportional to exp(S[a]) (uniform prior over the 20 amino acids), computed over the
  consensus columns the codon aligned to (N of them)
- a codon is called a when max P >= PMIN and N >= NMIN; otherwise it is reported as '?'
A codon call is non-standard (NS) when the called amino acid differs from the standard code. For TGA,
"standard" means not translated, so any confident amino-acid call for TGA is NS (UGA->Trp/Gly reassignments).

Clade level: S and N are summed over the genomes of each genus and each family, and the call is repeated on the
pooled evidence. The script reports NS clade calls with the supporting genome counts.

usage: gc_infer.py screen.npz genomes.txt taxmap.tsv meta_bac.tsv.gz meta_ar.tsv.gz OUTPREFIX
"""
import sys, numpy as np, pandas as pd
NPZ, GEN, TAX, MB, MA, OUT, HMMF = sys.argv[1:8]
AA = 'ACDEFGHIKLMNPQRSTVWY'
B = 'TCAG'; COD = [a + b + c for a in B for b in B for c in B]
STD = dict(zip(COD, 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
PMIN, NMIN = 0.9999, 5

d = np.load(NPZ); S, N, C = d['S'].astype(np.float64), d['N'], d['C']
genomes = open(GEN).read().split()
tax = pd.read_csv(TAX, sep='\t').set_index('accession').reindex(genomes)
meta = pd.concat([pd.read_csv(f, sep='\t', usecols=['accession', 'checkm2_completeness', 'checkm2_contamination',
                  'gc_percentage', 'genome_size', 'coding_density', 'ncbi_translation_table', 'ncbi_genome_category'],
                  low_memory=False) for f in (MB, MA)]).set_index('accession').reindex(genomes)

# Codetta normalisation: P(col k | a) = e_k(a) / sum_j e_j(a). The column population is (nearly) the same
# set of marker consensus columns in every genome, so sum_j e_j(a) is approximated by N_total * zbar[a],
# where zbar is the mean emission of amino acid a over the consensus columns of the marker HMMs.
import pyhmmer
_E = []
with pyhmmer.plan7.HMMFile(HMMF) as f:
    for h in f:
        em = np.array(h.match_emissions)[1:, :20]
        _E.append(em[em.max(1) >= 0.5])
_E = np.vstack(_E); LOGZ = np.log(_E.mean(0))

def call(Sv, Nv):
    """Sv: (..., 20) summed log emissions; Nv: (...) column counts.
    Returns best aa index (20 = nonspecific), posterior of best, confident flag."""
    L = Sv - Nv[..., None] * LOGZ[None if Sv.ndim == 2 else None]          # aa models minus common -m*log(Ntot)
    L = np.concatenate([L, np.zeros(L.shape[:-1] + (1,))], axis=-1)        # nonspecific model (reference = 0)
    m = L.max(-1, keepdims=True)
    p = np.exp(L - m); p /= p.sum(-1, keepdims=True)
    best = p.argmax(-1); pmax = p.max(-1)
    ok = (pmax >= PMIN) & (Nv >= NMIN) & (best < 20)
    return best, pmax, ok

best, pmax, ok = call(S, N)
std_idx = np.array([AA.find(STD[c]) for c in COD])  # -1 for stop
ns = ok & (best != std_idx[None, :])
# genome-level table of NS calls
rows = []
gi, ci = np.nonzero(ns)
for g, c in zip(gi, ci):
    s = S[g, c]; a = best[g, c]; sa = std_idx[c]
    llr = s[a] - (s[sa] if sa >= 0 else np.sort(s)[-2])
    rows.append(dict(genome=genomes[g], codon=COD[c], std=STD[COD[c]], call=AA[a], n=int(N[g, c]), pmax=float(pmax[g, c]),
                     llr=float(llr), n_cons_call=int(C[g, c, a]),
                     genus=tax.at[genomes[g], 'genus'], family=tax.at[genomes[g], 'family'], order=tax.at[genomes[g], 'order'],
                     phylum=tax.at[genomes[g], 'phylum']))
gt = pd.DataFrame(rows)
gt = gt.merge(meta, left_on='genome', right_index=True, how='left')
gt.to_csv(OUT + '_genome_ns.tsv', sep='\t', index=False)
print('genomes', len(genomes), 'NS genome-codon calls', len(gt), 'genomes with any NS', gt.genome.nunique())
# overall calling stats
print('called fraction of genome-codons (sense, N>=5):', ok[:, std_idx >= 0].mean().round(4))
acc = (best == std_idx[None, :])[ok & (std_idx[None, :] >= 0)].mean()
print('agreement with standard code among confident sense calls:', round(float(acc), 6))

# clade pooling
for lev in ['genus', 'family']:
    labs = tax[lev].fillna('NA').values
    df = []
    for lab in np.unique(labs):
        idx = np.where(labs == lab)[0]
        if len(idx) < 2: continue
        Sp = S[idx].sum(0); Np = N[idx].sum(0)
        b, p, o = call(Sp, Np)
        for c in range(64):
            if o[c] and b[c] != std_idx[c]:
                # how many member genomes individually support / contradict
                sup = int((ok[idx, c] & (best[idx, c] == b[c])).sum())
                con = int((ok[idx, c] & (best[idx, c] == std_idx[c])).sum())
                df.append(dict(level=lev, clade=lab, n_genomes=len(idx), codon=COD[c], std=STD[COD[c]], call=AA[b[c]],
                               n_cols=int(Np[c]), genomes_support=sup, genomes_standard=con,
                               phylum=tax.iloc[idx[0]]['phylum'], cls=tax.iloc[idx[0]]['class'], order=tax.iloc[idx[0]]['order']))
    pd.DataFrame(df).to_csv(OUT + f'_{lev}_ns.tsv', sep='\t', index=False)
    print(lev, 'NS clade-codon calls:', len(df))
