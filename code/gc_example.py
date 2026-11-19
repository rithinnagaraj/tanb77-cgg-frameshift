#!/usr/bin/env python
"""Alignment example: the most frequently hit conserved-Gly marker site. Shows CAG-508 genomes that carry CGG
at the site next to reference Clostridia, with a +/-8 residue window. Writes res/example_alignment.tsv."""
import pandas as pd, numpy as np, pyhmmer, os, random
s = pd.read_csv('res/cag508_cgg_sites_gc.tsv', sep='\t')
nat = s[(s.gene_gc <= 0.333) & (s.ref_fracG >= 0.95)]
top = nat.groupby(['marker', 'col']).size().sort_values(ascending=False).head(3)
print(top)
t = pd.read_csv(os.path.expanduser('~/work/taxmap.tsv'), sep='\t').set_index('accession')
root = 'data/bac120_marker_genes_reps_r226'
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
with pyhmmer.plan7.HMMFile('data/markers.hmm') as f:
    for h in f:
        a = h.accession or h.name; a = a.decode() if isinstance(a, bytes) else a; hm[a.split('.')[0]] = h
alpha = pyhmmer.easel.Alphabet.amino(); al = pyhmmer.plan7.TraceAligner()
rows = []
for (mk, col), cnt in top.items():
    P = rfa(f'{root}/faa/{mk}.faa'); N = rfa(f'{root}/fna/{mk}.fna')
    foc = list(nat[(nat.marker == mk) & (nat.col == col)].genome.unique())[:8]
    random.seed(1)
    refpool = [g for g in P if g in t.index and t.at[g, 'class'] == 'c__Clostridia' and t.at[g, 'order'] != 'o__TANB77']
    ref = random.sample(refpool, 6)
    others = [g for g in P if g in t.index and t.at[g, 'family'] == 'f__CAG-508' and g not in foc][:4]
    use = foc + others + ref
    seqs = [pyhmmer.easel.TextSequence(name=g.encode(), sequence=P[g].replace('*', '')).digitize(alpha) for g in use]
    blk = pyhmmer.easel.DigitalSequenceBlock(alpha, seqs)
    msa = al.align_traces(hm[mk.split('.')[0]], blk, al.compute_traces(hm[mk.split('.')[0]], blk), all_consensus_cols=True, trim=False)
    isx = [c == 'x' for c in msa.reference]
    colpos = [i for i, x in enumerate(isx) if x][col]
    for g, row in zip(use, msa.alignment):
        # residue index at colpos
        ri = sum(1 for c in row[:colpos] if c.isalpha())
        cod = N[g][3 * ri:3 * ri + 3] if row[colpos].isalpha() else '---'
        win = ''.join(c for i, c in enumerate(row[colpos - 12:colpos + 13]) if isx[colpos - 12 + i] or c.isalpha()).upper().replace('.', '')
        grp = 'CAG-508 (CGG)' if g in foc else ('CAG-508 (other)' if g in others else 'other Clostridia')
        rows.append(dict(marker=mk, col=col, group=grp, genome=g, genus=t.at[g, 'genus'], window=win, residue=row[colpos], codon=cod))
    print(mk, col, 'done')
df = pd.DataFrame(rows); df.to_csv('res/example_alignment.tsv', sep='\t', index=False)
print(df.to_string())
