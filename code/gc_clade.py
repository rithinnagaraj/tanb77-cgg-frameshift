#!/usr/bin/env python
"""Per-clade detail for one codon from the marker screen: per-genome evidence, pooled posteriors, and the
amino acid consensus counts of the columns the codon aligned to.
usage: gc_clade.py screen.npz genomes.txt taxmap.tsv markers.hmm CODON RANK NAME [RANK NAME ...]
"""
import sys, numpy as np, pandas as pd, pyhmmer
NPZ, GEN, TAX, HMMF, CODON = sys.argv[1:6]
pairs = list(zip(sys.argv[6::2], sys.argv[7::2]))
AA = 'ACDEFGHIKLMNPQRSTVWY'; B = 'TCAG'; COD = [a + b + c for a in B for b in B for c in B]
d = np.load(NPZ); S, N, C = d['S'].astype(float), d['N'], d['C']
genomes = open(GEN).read().split()
tax = pd.read_csv(TAX, sep='\t').set_index('accession').reindex(genomes)
E = []
with pyhmmer.plan7.HMMFile(HMMF) as f:
    for h in f:
        em = np.array(h.match_emissions)[1:, :20]; E.append(em[em.max(1) >= 0.5])
LOGZ = np.log(np.vstack(E).mean(0))
ci = COD.index(CODON)
def post(Sv, n):
    L = np.append(Sv - n * LOGZ, 0.0); p = np.exp(L - L.max()); return p / p.sum()
for rank, name in pairs:
    idx = np.where(tax[rank].values == name)[0]
    if len(idx) == 0: print(rank, name, 'none'); continue
    Sp = S[idx, ci].sum(0); Np = N[idx, ci].sum(); Cp = C[idx, ci].sum(0)
    p = post(Sp, Np); top = np.argsort(-p)[:4]
    nz = (N[idx, ci] > 0).sum()
    # per-genome posterior of Gly vs Arg among genomes with >=1 column
    g_gly = sum(1 for i in idx if N[i, ci] > 0 and post(S[i, ci], N[i, ci])[AA.index('G')] > 0.5)
    g_arg = sum(1 for i in idx if N[i, ci] > 0 and post(S[i, ci], N[i, ci])[AA.index('R')] > 0.5)
    cons = ', '.join(f'{AA[j]}:{Cp[j]}' for j in np.argsort(-Cp)[:6] if Cp[j] > 0)
    # codon usage of the codon in marker consensus columns relative to all codons
    tot = N[idx].sum()
    print(f'{rank} {name}: genomes={len(idx)} with_cols={nz} cols={Np} usage_per_10k={1e4*Np/max(tot,1):.2f} '
          f'pooled={[(AA+"*")[j]+":"+format(p[j],".3g") for j in top]} genomes_Gly>0.5={g_gly} genomes_Arg>0.5={g_arg} consensus=[{cons}]')
