#!/usr/bin/env python
"""Figures for the executive summary (static PNG, light theme).
Palette: categorical slot 1 (blue #2a78d6) for the focal clade; neutral grays for controls; text in ink tokens."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
BLUE = '#2a78d6'; ORANGE = '#eb6834'; AQUA = '#1baf7a'; GRAY = '#b9b7b0'; GRAY2 = '#8d8b85'
INK = '#0b0b0b'; INK2 = '#52514e'; SURF = '#fcfcfb'; GRID = '#e6e5e0'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.edgecolor': GRAY2, 'axes.labelcolor': INK2,
                     'xtick.color': INK2, 'ytick.color': INK2, 'axes.titlecolor': INK, 'axes.titlesize': 10,
                     'axes.titleweight': 'bold', 'figure.facecolor': 'white', 'axes.facecolor': SURF,
                     'axes.spines.top': False, 'axes.spines.right': False})
OUT = 'figs/'
import os; os.makedirs(OUT, exist_ok=True)

def bars(ax, labels, vals, colors, fmt='{:.1f}%', ymax=None):
    x = np.arange(len(vals))
    ax.bar(x, vals, width=0.6, color=colors, edgecolor=SURF, linewidth=2)
    for i, v in enumerate(vals):
        ax.text(i, v + (ymax or max(vals)) * 0.02, fmt.format(v), ha='center', va='bottom', color=INK, fontsize=8.5)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8)
    ax.yaxis.grid(True, color=GRID, lw=0.8); ax.set_axisbelow(True)
    if ymax: ax.set_ylim(0, ymax)

# ---------- Fig 1: genome-level frame test + next nucleotide
f = pd.read_csv('res/frames_genome.tsv', sep='\t'); f['delta'] = f.score_f1 - f.score_f0
g = f.groupby(['clade', 'kind']).agg(n=('delta', 'size'), better=('delta', lambda x: 100 * (x > 0).mean()))
order = [('CAG-508', 'CGG'), ('CAG-508', 'GGN_ctrl'), ('other Clostridia', 'CGG'), ('other Clostridia', 'GGN_ctrl')]
lab = [f'{a}\n{b}\nn={int(g.loc[o, "n"]):,}' for (a, b), o in zip([('CAG-508', 'CGG'), ('CAG-508', 'GGN ctrl'), ('Other', 'CGG'), ('Other', 'GGN ctrl')], order)]
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.9), gridspec_kw={'width_ratios': [1.2, 1]})
bars(axs[0], lab, [g.loc[o, 'better'] for o in order], [BLUE, GRAY, GRAY, GRAY], ymax=100)
axs[0].tick_params(axis='x', labelsize=7)
axs[0].set_ylabel('% sites where the +1 frame\nmatches the homolog better')
axs[0].set_title('A  +1 frame restores homology', loc='left', fontsize=9)
c = f[f.kind == 'CGG'].groupby('clade').next_nt.value_counts(normalize=True).unstack().fillna(0) * 100
nts = ['A', 'T', 'G', 'C']; x = np.arange(4); w = 0.36
axs[1].bar(x - w / 2, [c.loc['CAG-508', n] for n in nts], w, color=BLUE, edgecolor=SURF, lw=2, label='CAG-508')
axs[1].bar(x + w / 2, [c.loc['other Clostridia', n] for n in nts], w, color=GRAY, edgecolor=SURF, lw=2, label='Other Clostridia')
axs[1].set_xticks(x); axs[1].set_xticklabels(['CGG·A', 'CGG·U', 'CGG·G', 'CGG·C'])
axs[1].set_ylabel('% of CGG sites'); axs[1].yaxis.grid(True, color=GRID); axs[1].set_axisbelow(True)
axs[1].legend(frameon=False, fontsize=8); axs[1].set_title('B  Base following CGG', loc='left', fontsize=9)
plt.tight_layout(); plt.savefig(OUT + 'fig1_frame_test.png', dpi=220); plt.close()

# ---------- Fig 2: alignment example (text figure)
ex = [
 ('TIGR02386 (RNA polymerase β′, col 696)', [
   ('Reference Clostridia (6/6)', 'NPIYMMAQSGARG NMNQLRQIAGMR', 'GGT'),
   ('CAG-508, annotated frame', 'NPIYMMAQSGARR KHEPI* ', 'CGG'),
   ('CAG-508, read +1 after C', 'NPIYMMAQSGARG NMNQLRQIAGMR', 'C|GGA'),
   ('CAG-508 nucleotides', '...TCA GGA GCA AGA [C]GGA AAC ATG AAC CAG TTA AGA CAA ATT GCA GGT', '')]),
 ('TIGR00580 (Mfd, col 808; 84 genomes, 13 genera)', [
   ('Reference Clostridia', 'RDLEIRGAGNLLG AEQHGHMDKVGY', 'GGC'),
   ('CAG-508, annotated frame', 'RDLEIRGAGSLLR RNTTWSYGTSWI*', 'CGG'),
   ('CAG-508, read +1 after C', 'RDLEIRGAGSLLG EIQHGHMEQVGY', 'C|GGA'),
   ('CAG-508 nucleotides', '...GGA AGT TTA CTT [C]GGA GAA ATA CAA CAT GGT CAT ATG GAA CAA GTT', '')])]
fig, ax = plt.subplots(figsize=(7.2, 2.9)); ax.axis('off'); y = 0.95
for title, rows in ex:
    ax.text(0, y, title, fontsize=9, fontweight='bold', color=INK, transform=ax.transAxes); y -= 0.085
    for name, seq, cod in rows:
        col = BLUE if 'read +1' in name else (ORANGE if 'annotated' in name else INK)
        ax.text(0.0, y, name, fontsize=7.8, color=INK2, transform=ax.transAxes)
        ax.text(0.29, y, seq, fontsize=7.8, family='DejaVu Sans Mono', color=col, transform=ax.transAxes)
        y -= 0.075
    y -= 0.04
ax.text(0, 0.0, 'Orange: Prodigal/GTDB annotation (CGG read in frame, homology lost, early stop). Blue: the same DNA read after a +1 shift at the extra C.',
        fontsize=7, color=INK2, transform=ax.transAxes)
plt.savefig(OUT + 'fig2_alignment.png', dpi=220, bbox_inches='tight'); plt.close()

# ---------- Fig 3: tRNA loss and frameshift signature across TANB77 families + control
gw = pd.read_csv('res/genomewide_annot.tsv', sep='\t')
t = pd.read_csv(os.path.expanduser('~/work/taxmap.tsv'), sep='\t').set_index('accession')
pr = pd.read_csv('res/profile_annot.tsv', sep='\t').set_index('genome')
gw['family'] = gw.genome.map(t.family); gw['reader'] = gw.genome.map(pr.has_CCG_tRNA | pr.has_UCG_tRNA)
gw['fam2'] = np.where(gw.grp == 'other Clostridia', 'Other Clostridia', gw.family.str[3:])
fams = ['RGIG9151', 'CAG-465', 'UBA1234', 'CAG-508', 'Other Clostridia']
s = gw.groupby('fam2').agg(n=('genome', 'size'), reader=('reader', 'mean'), term=('CGG_frac_term', 'median'), shift=('CGG_shiftfrac', 'median')).reindex(fams)
fig, axs = plt.subplots(1, 3, figsize=(7.6, 2.9))
cols = [GRAY2, BLUE, BLUE, BLUE, GRAY]
SHORT = {'RGIG9151': 'RGIG\n9151', 'CAG-465': 'CAG-\n465', 'UBA1234': 'UBA\n1234', 'CAG-508': 'CAG-\n508', 'Other Clostridia': 'Other\nClostr.'}
labs = [SHORT[f] + '\n(' + str(int(s.loc[f, 'n'])) + ')' for f in fams]
bars(axs[0], labs, list(100 * s.reader), cols, fmt='{:.0f}%', ymax=110); axs[0].set_title('A  Has a CGG-reading tRNA\n(anticodon CCG or UCG)', loc='left', fontsize=8.5)
bars(axs[1], labs, list(100 * s.term), cols, fmt='{:.0f}%', ymax=55); axs[1].set_title('B  CGGs in last 15 codons\nof an annotated gene', loc='left', fontsize=8.5)
bars(axs[2], labs, list(100 * s['shift']), cols, fmt='{:.0f}%', ymax=85); axs[2].set_title('C  …with a +1 ORF continuing\n≥60 codons past the stop', loc='left', fontsize=8.5)
for a in axs: a.tick_params(axis='x', labelsize=6.2, rotation=0)
plt.tight_layout(); plt.savefig(OUT + 'fig3_trna_link.png', dpi=220); plt.close()

# ---------- Fig 4: codon specificity (rate per source codon, log scale) and rare-CGG null
AA = 'ACDEFGHIKLMNPQRSTVWY'; COD = [a + b + c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG']
d = np.load('screen_all.npz'); N = d['N']; C = d['C']
gg = open('screen_all.npz.genomes.txt').read().split(); tt = t.reindex(gg)
iG = AA.index('G')
groups = {'CAG-508': tt.family.values == 'f__CAG-508', 'Oscillospirales': tt.order.values == 'o__Oscillospirales',
          'Lachnospirales': tt.order.values == 'o__Lachnospirales', 'Christensenellales': tt.order.values == 'o__Christensenellales'}
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.8), gridspec_kw={'width_ratios': [1, 1.25]})
xs = np.arange(4); w = 0.2
for j, (name, m) in enumerate(groups.items()):
    vals = []
    for x3, src in [('G', 'GGG'), ('A', 'GGA'), ('T', 'GGT'), ('C', 'GGC')]:
        k = C[m, COD.index('CG' + x3), iG].sum(); n = N[m, COD.index(src)].sum(); vals.append(max(k / n, 1e-6))
    axs[0].bar(xs + (j - 1.5) * w, vals, w, color=BLUE if name == 'CAG-508' else [GRAY, GRAY2, '#d6d4ce'][j - 1], edgecolor=SURF, lw=1, label=name)
axs[0].set_yscale('log'); axs[0].set_xticks(xs); axs[0].set_xticklabels(['CGG\nper GGG', 'CGA\nper GGA', 'CGU\nper GGU', 'CGC\nper GGC'], fontsize=7.5)
axs[0].set_ylabel('events per source codon'); axs[0].legend(frameon=False, fontsize=6.3, loc='upper right'); axs[0].set_ylim(1e-6, 0.3); axs[0].yaxis.grid(True, color=GRID); axs[0].set_axisbelow(True)
axs[0].set_title('A  Codon-specific (439×), not\nsequencing error', loc='left', fontsize=8.5)
rows = []
fam = tt.family.fillna('NA').values; ci = COD.index('CGG')
nl = pd.read_csv('res/null_rare_codons.tsv', sep='\t'); nl = nl[nl.codon == 'CGG']
for f_ in nl.family:
    m = fam == f_; cc = C[m, ci].sum(0); tot = cc.sum()
    rows.append((f_[3:], 100 * cc[AA.index('R')] / tot, 100 * cc[iG] / tot, 100 * cc[AA.index('W')] / tot, int(N[m, ci].sum())))
r = pd.DataFrame(rows, columns=['family', 'Arg', 'Gly', 'Trp', 'n']).sort_values('Gly')
r['other'] = 100 - r.Arg - r.Gly - r.Trp
yy = np.arange(len(r)); left = np.zeros(len(r))
for k, cl in [('Arg', GRAY), ('Gly', BLUE), ('Trp', AQUA), ('other', '#e3e1db')]:
    axs[1].barh(yy, r[k], left=left, color=cl, edgecolor=SURF, lw=1, label=k, height=0.7); left += r[k].values
axs[1].set_yticks(yy); axs[1].set_yticklabels([f'{a} ({b})' for a, b in zip(r.family, r.n)], fontsize=6.5)
axs[1].set_xlabel('% of CGG columns by consensus residue'); axs[1].legend(frameon=False, fontsize=6.5, ncol=4, loc='lower center', bbox_to_anchor=(0.45, 1.0))
axs[1].set_title('B  All 13 families with rare CGG (<5/10k)', loc='left', fontsize=8.5, pad=16)
plt.tight_layout(); plt.savefig(OUT + 'fig4_specificity.png', dpi=220); plt.close()
# ---------- Fig 5: long-read replication + tRNA-Gly variant
lr = pd.read_csv('res/genomewide_lr.tsv', sep='\t')
lr['CGG'] = 100 * lr.CGG_shift / lr.CGG_term.clip(lower=1); lr['AGG'] = 100 * lr.AGG_shift / lr.AGG_term.clip(lower=1); lr['GGG'] = 100 * lr.GGG_shift / lr.GGG_term.clip(lower=1)
lr['term'] = 100 * lr.CGG_term / lr.CGG_n
gv = pd.read_csv('res/glyvar_annot.tsv', sep='\t')
gv['fam'] = np.where(gv.grp == 'other Clostridia', 'Other Clostridia', gv.family.str[3:])
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.7), gridspec_kw={'width_ratios': [1, 1.15]})
rng = np.random.default_rng(1)
for i, (k, cl) in enumerate([('CGG', BLUE), ('AGG', GRAY2), ('GGG', GRAY)]):
    xs_ = i + rng.uniform(-0.18, 0.18, len(lr))
    axs[0].scatter(xs_, lr[k], s=16, color=cl, edgecolor=SURF, linewidth=0.6, zorder=3)
    axs[0].plot([i - 0.28, i + 0.28], [lr[k].median()] * 2, color=INK, lw=1.5, zorder=4)
axs[0].set_xticks([0, 1, 2]); axs[0].set_xticklabels(['CGG', 'AGG\n(control)', 'GGG\n(control)'])
axs[0].set_ylabel('% of gene-tail sites with a\n+1 ORF continuing ≥60 codons'); axs[0].set_ylim(0, 100)
axs[0].yaxis.grid(True, color=GRID); axs[0].set_axisbelow(True)
axs[0].set_title('A  23 complete CAG-508 genomes\n(19 PacBio HiFi, 6 projects)', loc='left', fontsize=8.5)
fams = ['CAG-508', 'UBA1234', 'CAG-465', 'RGIG9151', 'Other Clostridia']
sv = gv.groupby('fam').agg(n=('genome', 'size'), var=('has_var', 'mean'), rd=('n_CGG_reader', lambda x: (x > 0).mean())).reindex(fams)
yy = np.arange(len(fams))[::-1]; h = 0.36
axs[1].barh(yy + h / 2, 100 * sv['var'], h, color=BLUE, edgecolor=SURF, lw=1, label='tRNA-Gly(UCC) G37 variant')
axs[1].barh(yy - h / 2, 100 * sv['rd'], h, color=GRAY, edgecolor=SURF, lw=1, label='CGG-reading tRNA (CCG/UCG)')
axs[1].set_yticks(yy); axs[1].set_yticklabels([f'{f} ({int(sv.loc[f, "n"])})' for f in fams], fontsize=7.5)
axs[1].set_xlabel('% of genomes'); axs[1].set_xlim(0, 105); axs[1].xaxis.grid(True, color=GRID); axs[1].set_axisbelow(True)
axs[1].legend(frameon=False, fontsize=6.8, loc='lower center', bbox_to_anchor=(0.4, 1.0), ncol=1)
axs[1].set_title('B  tRNA swap', loc='left', fontsize=8.5, pad=30)
plt.tight_layout(); plt.savefig(OUT + 'fig5_longread_trna.png', dpi=220); plt.close()

# ---------- Fig 0: model schematic
fig, ax = plt.subplots(figsize=(7.2, 1.9)); ax.axis('off'); ax.set_xlim(0, 100); ax.set_ylim(0, 30)
nts = list('GCAAGACGGAAACATGAAC'); x0 = 6; dx = 4.2
for i, ch in enumerate(nts):
    c = BLUE if 6 <= i <= 9 else INK
    ax.text(x0 + i * dx, 18, ch, family='DejaVu Sans Mono', fontsize=12, ha='center', color=c, fontweight='bold' if 6 <= i <= 9 else 'normal')
ax.text(1, 18, "5'", fontsize=9, color=INK2, va='center'); ax.text(x0 + len(nts) * dx, 18, "3'", fontsize=9, color=INK2)
# annotated frame (orange) codons
for j, (a, lab) in enumerate([(0, 'A'), (3, 'R'), (6, 'R'), (9, 'K'), (12, 'H'), (15, 'E')]):
    ax.add_patch(FancyBboxPatch((x0 + a * dx - 1.8, 23.2), 3 * dx - 0.6, 3.2, boxstyle='round,pad=0.1', fc='none', ec=ORANGE, lw=1))
    ax.text(x0 + (a + 1) * dx, 24.8, lab, fontsize=8, ha='center', va='center', color=ORANGE)
ax.text(x0 + 19 * dx + 1, 24.8, 'annotated frame (junk → stop)', fontsize=7.5, color=ORANGE, va='center')
for j, (a, w, lab) in enumerate([(0, 3, 'A'), (3, 3, 'R'), (6, 4, 'G'), (10, 3, 'N'), (13, 3, 'M'), (16, 3, 'N')]):
    ax.add_patch(FancyBboxPatch((x0 + a * dx - 1.8, 10.2), w * dx - 0.6, 3.2, boxstyle='round,pad=0.1', fc='none', ec=BLUE, lw=1.2 if w == 4 else 0.8))
    ax.text(x0 + (a + w / 2 - 0.5) * dx, 11.8, lab, fontsize=8, ha='center', va='center', color=BLUE)
ax.text(x0 + 19 * dx + 1, 11.8, 'proposed: C-GGA read as one Gly (+1)', fontsize=7.5, color=BLUE, va='center')
ax.text(x0 + 7.5 * dx, 3.5, "tRNA-Gly variant anticodon loop  …C-U-[U34 C35 C36]-G37-A…  →  G37 can pair the C of C-GGA",
        fontsize=7.5, color=INK2, ha='center')
plt.savefig(OUT + 'fig0_model.png', dpi=220, bbox_inches='tight'); plt.close()
print('figures done')
