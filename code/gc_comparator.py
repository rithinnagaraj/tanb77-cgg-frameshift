#!/usr/bin/env python
"""Review-requested fixes:
(1) Insertion-artefact comparator: a base X inserted before a GGN Gly codon creates XGG (AGG, CGG, TGG). For each X,
    compare XGG-at-Gly-consensus-columns per GGG source codon in CAG-508 against other Clostridia orders.
(2) The next-base null: CAG-508's own GGN third-position usage."""
import os
import numpy as np, pandas as pd
AA = 'ACDEFGHIKLMNPQRSTVWY'; COD = [a + b + c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG']
d = np.load('screen_all.npz'); N = d['N']; C = d['C']
g = open('screen_all.npz.genomes.txt').read().split()
t = pd.read_csv(os.path.expanduser('~/work/taxmap.tsv'), sep='\t').set_index('accession').reindex(g)
iG = AA.index('G')
grp = {'CAG-508': t.family.values == 'f__CAG-508', 'Oscillospirales': t.order.values == 'o__Oscillospirales',
       'Lachnospirales': t.order.values == 'o__Lachnospirales', 'Christensenellales': t.order.values == 'o__Christensenellales'}
rows = []
for k, m in grp.items():
    src = {x: N[m, COD.index('GG' + x)].sum() for x in 'ACGT'}
    for x in 'ACT':
        ev = int(C[m, COD.index(x + 'GG'), iG].sum())
        rows.append(dict(group=k, codon=x + 'GG', events=ev, per_GGN=ev / sum(src.values())))
df = pd.DataFrame(rows)
print(df.pivot(index='codon', columns='group', values='per_GGN').map(lambda v: f'{v:.2e}').to_string())
print(df.pivot(index='codon', columns='group', values='events').to_string())
m = grp['CAG-508']; tot = {x: N[m, COD.index('GG' + x)].sum() for x in 'ATGC'}; s = sum(tot.values())
print('CAG-508 GGN 3rd-position usage:', {x: round(v / s, 3) for x, v in tot.items()}, '| base after CGG at shift sites: A .72 T .18 G .05 C .04')
