#!/usr/bin/env python
"""Assemble per-genus figure data for TANB77, ordered by the GTDB bac120 r226 tree, plus an alignment example."""
import os
import pandas as pd, numpy as np, re, sys
from ete3 import Tree
t = pd.read_csv(os.path.expanduser('~/work/taxmap.tsv'), sep='\t')
tan = t[t.order == 'o__TANB77']
tree = Tree('data/bac120_r226.tree', format=1, quoted_node_names=True)
keep = set(tan.accession)
leaves = [l for l in tree.iter_leaves() if l.name in keep]
anc = tree.get_common_ancestor(leaves)
order = [l.name for l in anc.iter_leaves() if l.name in keep]
gen_order = list(dict.fromkeys(tan.set_index('accession').loc[order, 'genus']))
prof = pd.read_csv('res/profile_annot.tsv', sep='\t')
prof = prof[prof.genome.isin(keep)].merge(tan[['accession', 'genus', 'family']], left_on='genome', right_on='accession', suffixes=('', '_t'))
g1 = prof.groupby('genus_t').agg(n=('genome', 'size'), cgg10k=('CGG_per_10k', 'median'), ccg_trna=('has_CCG_tRNA', 'mean'))
s = pd.read_csv('res/cag508_cgg_sites_gc.tsv', sep='\t'); nat = s[s.gene_gc <= 0.333]
g2 = nat.groupby('genus').agg(sitesG=('ref_top', lambda x: (x == 'G').sum()), sitesR=('ref_top', lambda x: (x == 'R').sum()))
mel = pd.read_csv('lit/melnykov/bacteria_summary.csv.gz')
COD = [a + b + c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG']
mel['CGG'] = mel['inferred genetic code'].str[COD.index('CGG')]; mel['nov'] = mel.assembly.str.split('.').str[0]
md = pd.read_csv(os.path.expanduser('~/data/gtdb/bac120_metadata_r226.tsv.gz'), sep='\t', usecols=['accession', 'gtdb_taxonomy'], low_memory=False)
md = md[md.gtdb_taxonomy.str.contains('o__TANB77')]
md['nov'] = md.accession.str[3:].str.split('.').str[0]; md['genus'] = md.gtdb_taxonomy.str.split(';').str[5]
mj = mel.merge(md, on='nov')
g3 = mj.groupby('genus').CGG.value_counts().unstack(fill_value=0)
fam = tan.groupby('genus').family.first()
df = pd.DataFrame(index=gen_order).join(fam).join(g1).join(g2).join(g3.add_prefix('kaci_')).fillna(0)
df.index.name = 'genus'
df.to_csv('res/fig_tanb77_genera.tsv', sep='\t')
print(df.head(60).to_string())
