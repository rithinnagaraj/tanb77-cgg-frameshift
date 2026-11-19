#!/usr/bin/env python
"""Statistical checks requested by the adversarial review (run on the marker-screen arrays).

1. Calibrated null: for every GTDB family and every *rare* sense codon (usage < 5 per 10k marker consensus
   columns, at least 30 columns pooled), compute the fraction of the codon's columns whose consensus is the
   single most frequent amino acid that is NOT the standard one. This gives the empirical tail probability of
   CAG-508 CGG -> 0.85 Gly.
2. Rate model: per clade, count codon-at-Gly-column events per source codon.
   CGG@G / GGG  versus  CGA@G / GGA, CGT@G / GGT, CGC@G / GGC  (as Poisson rates, with exact CIs).
3. Gly vs Ala vs Ser posterior margins for pooled CAG-508 CGG, and the RGIG3102 control.
4. Pseudo-replication:
   (a) one genome per genus (best-evidence genome excluded; a random genome chosen),
   (b) leave-out CAG-269,
   (c) the number of distinct (marker, column) sites (from the site table).
"""
import os
import numpy as np, pandas as pd, pyhmmer, sys
from scipy import stats
AA = 'ACDEFGHIKLMNPQRSTVWY'; B = 'TCAG'; COD = [a + b + c for a in B for b in B for c in B]
STD = dict(zip(COD, 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
d = np.load('screen_all.npz'); S = d['S'].astype(float); N = d['N']; C = d['C']
g = open('screen_all.npz.genomes.txt').read().split()
t = pd.read_csv(os.path.expanduser('~/work/taxmap.tsv'), sep='\t').set_index('accession').reindex(g)
E = []
with pyhmmer.plan7.HMMFile('data/markers.hmm') as f:
    for h in f:
        em = np.array(h.match_emissions)[1:, :20]; E.append(em[em.max(1) >= 0.5])
LOGZ = np.log(np.vstack(E).mean(0))
def post(Sv, n):
    L = np.append(Sv - n * LOGZ, 0.0); p = np.exp(L - L.max()); return p / p.sum()
out = open('res/checks.txt', 'w')
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')

# ---- 1. calibrated null over rare codons in all families
fam = t.family.fillna('NA').values
codtot = N.sum(1)
rows = []
for f_ in np.unique(fam):
    idx = np.where(fam == f_)[0]
    if len(idx) < 5: continue
    tot = N[idx].sum()
    for ci, cod in enumerate(COD):
        if STD[cod] == '*': continue
        n = N[idx, ci].sum()
        if n < 30: continue
        usage = 1e4 * n / tot
        if usage >= 5: continue
        cc = C[idx, ci].sum(0); std_i = AA.index(STD[cod])
        alt = cc.copy(); alt[std_i] = -1; ai = int(alt.argmax())
        rows.append(dict(family=f_, codon=cod, std=STD[cod], n=int(n), usage=usage, frac_std=cc[std_i] / cc.sum(),
                         top_alt=AA[ai], frac_alt=cc[ai] / cc.sum()))
nl = pd.DataFrame(rows).sort_values('frac_alt', ascending=False)
nl.to_csv('res/null_rare_codons.tsv', sep='\t', index=False)
P('== 1. Calibrated null: rare (usage<5/10k, >=30 cols) family x codon combos:', len(nl))
P(nl.head(12).to_string())
x = nl[(nl.family != 'f__CAG-508')]
P('combos (excluding CAG-508) with frac_alt >= 0.5:', int((x.frac_alt >= 0.5).sum()), 'of', len(x),
  '; max frac_alt =', round(float(x.frac_alt.max()), 3))
P('empirical P(frac_alt >= CAG-508 value) =', round(float((x.frac_alt >= float(nl[(nl.family == 'f__CAG-508') & (nl.codon == 'CGG')].frac_alt.iloc[0])).mean()), 5)
  if len(nl[(nl.family == 'f__CAG-508') & (nl.codon == 'CGG')]) else 'NA')

# ---- 2. rate model
P('\n== 2. Rate of X-at-Gly-consensus columns per source codon (events / source-codon columns)')
grp = {'CAG-508': t.family.values == 'f__CAG-508', 'other TANB77': (t.order.values == 'o__TANB77') & (t.family.values != 'f__CAG-508'),
       'Oscillospirales': t.order.values == 'o__Oscillospirales', 'Lachnospirales': t.order.values == 'o__Lachnospirales',
       'Christensenellales': t.order.values == 'o__Christensenellales'}
iG = AA.index('G')
for name, m in grp.items():
    parts = []
    for x3, src in [('G', 'GGG'), ('A', 'GGA'), ('T', 'GGT'), ('C', 'GGC')]:
        k = int(C[m, COD.index('CG' + x3), iG].sum()); nsrc = int(N[m, COD.index(src)].sum())
        lo, hi = stats.poisson.interval(0.95, k) if k > 0 else (0, 3.0)
        parts.append(f'CG{x3}@G/{src}={k}/{nsrc}={k / max(nsrc, 1):.2e}')
    P(name, ' | '.join(parts))
m = grp['CAG-508']
k1 = int(C[m, COD.index('CGG'), iG].sum()); n1 = int(N[m, COD.index('GGG')].sum())
k0 = sum(int(C[m, COD.index('CG' + x), iG].sum()) for x in 'ATC'); n0 = sum(int(N[m, COD.index('GG' + x)].sum()) for x in 'ATC')
rr = (k1 / n1) / (k0 / n0)
p = stats.binomtest(k1, k1 + k0, n1 / (n1 + n0), alternative='greater').pvalue
P(f'CAG-508 rate ratio CGG@G per GGG vs CGH@G per GGH: {rr:.1f}x (binomial p={p:.2e}); k1={k1}, n1={n1}, k0={k0}, n0={n0}')

# ---- 3. margins and the RGIG3102 control
P('\n== 3. Posterior margins (pooled)')
for lab, m in [('f__CAG-508', t.family.values == 'f__CAG-508'), ('g__RGIG3102', t.genus.values == 'g__RGIG3102'),
               ('g__Peptacetobacter', t.genus.values == 'g__Peptacetobacter'), ('g__Onthovivens', t.genus.values == 'g__Onthovivens')]:
    ci = COD.index('CGG'); pp = post(S[m, ci].sum(0), N[m, ci].sum())
    lab2 = AA + '*'; top = np.argsort(-pp)[:4]
    P(lab, 'genomes', int(m.sum()), 'cols', int(N[m, ci].sum()), [(lab2[j], f'{pp[j]:.3g}') for j in top],
      'logL(G)-logL(A)=', round(float((S[m, ci].sum(0)[iG] - N[m, ci].sum() * LOGZ[iG]) - (S[m, ci].sum(0)[AA.index('A')] - N[m, ci].sum() * LOGZ[AA.index('A')])), 1))

# ---- 4. pseudo-replication
P('\n== 4. Pseudo-replication')
rng = np.random.default_rng(0); ci = COD.index('CGG')
cag = np.where(t.family.values == 'f__CAG-508')[0]
for rep in range(5):
    pick = [rng.choice(np.where(t.genus.values == gg)[0][np.isin(np.where(t.genus.values == gg)[0], cag)]) for gg in np.unique(t.genus.values[cag])]
    pick = np.array(pick); pp = post(S[pick, ci].sum(0), N[pick, ci].sum()); cc = C[pick, ci].sum(0)
    P(f'one-genome-per-genus rep{rep}: genomes={len(pick)} cols={int(N[pick, ci].sum())} P(G)={pp[iG]:.4f} Gcons={cc[iG]} Rcons={cc[AA.index("R")]}')
lo = cag[t.genus.values[cag] != 'g__CAG-269']
pp = post(S[lo, ci].sum(0), N[lo, ci].sum()); cc = C[lo, ci].sum(0)
P(f'leave-out CAG-269: genomes={len(lo)} cols={int(N[lo, ci].sum())} P(G)={pp[iG]:.6f} Gcons={cc[iG]} Rcons={cc[AA.index("R")]}')
s = pd.read_csv('res/cag508_cgg_sites_gc.tsv', sep='\t'); nat = s[s.gene_gc <= 0.333]
P('distinct (marker,column) sites, native copies:', nat.groupby(['marker', 'col']).ngroups, 'from', len(nat), 'occurrences;',
  'distinct sites with ref-majority Gly:', nat[nat.ref_top == 'G'].groupby(['marker', 'col']).ngroups,
  'Arg:', nat[nat.ref_top == 'R'].groupby(['marker', 'col']).ngroups)
P('native sites, ref-majority residue (top 6):', nat.ref_top.value_counts().head(6).to_dict())
out.close()
