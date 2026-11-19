#!/usr/bin/env python
"""Blind positive-control evaluation (sensitivity) and the per-job covariation distribution.

Ground truth comes from Rfam cmscan hits (GA thresholds) of the expected riboswitch family on the
sequences of the positive-control jobs.
- A job is detectable if at least MINTRUE of its sequences carry the expected riboswitch.
- A job counts as recovered if some motif has at least 50% of its instances overlapping (>= 30% of the
  Rfam hit) the expected riboswitch, and that motif passes the significance rule tested.
- A motif is "known" if more than 30% of its instances overlap any Rfam hit (this is used for the targets
  and the nulls).

usage: eval_controls.py JOBDIR_MANIFEST MOTIFS_TSV MOTIFS_FA JOBS_RFAM_TBL MOTIF_RFAM_TBL FAMILIES OUTPREFIX
"""
import sys, re, collections
import pandas as pd, numpy as np

MAN, MT, MFA, JTBL, MTBL, FAM, OUT = sys.argv[1:8]
EXPECT = {
 'POS_TPP': r'^TPP$', 'POS_FMN': r'^FMN$', 'POS_COBALAMIN': r'Cobalamin|AdoCbl', 'POS_NAD': r'^nadA$',
 'POS_MOCO': r'Molybdenum|Moco', 'POS_SAM': r'^SAM', 'POS_LYS': r'^Lysine$', 'POS_GLY': r'^Glycine$|glycine-GGAnGA',
 'POS_FLUORIDE': r'^Fluoride$', 'POS_PREQ1': r'PreQ1|preQ1', 'POS_MG': r'ykoK|M-box', 'POS_GLMS': r'^glmS$',
 'POS_PURINE': r'^Purine$|2dG', 'POS_LI': r'nhaA', 'POS_NICO': r'NiCo', 'POS_GUANIDINE': r'Guanidine',
 'POS_ZTP': r'ZMP-ZTP', 'POS_SAH': r'SAH', 'POS_HMPPP': r'HMP|thiS'}
MINTRUE = 3

def read_tbl(fn):
    """cmscan --fmt 2 tblout -> dict seqname -> list of (family, start, end)"""
    d = collections.defaultdict(list)
    for l in open(fn):
        if l.startswith('#'): continue
        p = l.split()
        if len(p) < 20: continue
        fam, seq = p[1], p[3]
        a, b = int(p[9]), int(p[10])
        strand, inc, olp = p[11], p[18], p[19]
        if inc != '!' or strand != '+' or olp == '=': continue
        d[seq].append((fam, min(a, b), max(a, b)))
    return d

fam = pd.read_csv(FAM, sep='\t'); grp = dict(zip(fam.ko, fam.group))
man = pd.read_csv(MAN, sep='\t')
mt = pd.read_csv(MT, sep='\t')
jt = read_tbl(JTBL)
mtb = read_tbl(MTBL)

# motif instances: header job|motif|uid|pos|len
inst = collections.defaultdict(list)
for l in open(MFA):
    if l[0] == '>':
        p = l[1:].strip().split('|')
        job, mo, uid, pos, ln = p[0], p[1], '|'.join(p[2:-2]), p[-2], p[-1]
        inst[(job, mo)].append((uid, int(pos), int(ln)))

def ov(a1, a2, b1, b2):
    return max(0, min(a2, b2) - max(a1, b1) + 1)

rows = []
for (job, mo), ins in inst.items():
    # overlap with any Rfam hit in the motif's own instances (motif-level novelty scan, coordinates 1..len)
    known = collections.Counter()
    for uid, pos, ln in ins:
        name = f'{job}|{mo}|{uid}|{pos}|{ln}'
        for f, a, b in mtb.get(name, []):
            if ov(1, ln, a, b) >= 0.3 * min(ln, b - a + 1): known[f] += 1
    top = known.most_common(1)
    rows.append(dict(job=job, motif=mo, known_fam=top[0][0] if top else '',
                     known_frac=(top[0][1] / len(ins)) if top else 0.0))
kn = pd.DataFrame(rows)
mt = mt.merge(kn, on=['job', 'motif'], how='left')
mt['group'] = mt.ko.map(lambda k: grp.get(k, 'NULL'))
mt['is_known'] = mt.known_frac > 0.3
mt.to_csv(OUT + '_motifs_annot.tsv', sep='\t', index=False)

# sensitivity on positive-control jobs
pos = man[man.ko.map(lambda k: grp.get(k, 'NULL')).str.startswith('POS')].copy()
res = []
for r in pos.itertuples():
    g = grp[r.ko]; pat = re.compile(EXPECT.get(g, '^$'))
    # job sequences carrying the expected riboswitch (job fasta coordinates)
    seqs = [l[1:].strip() for l in open(f'{sys.argv[1].rsplit("/",1)[0]}/{r.job}.fa') if l[0] == '>']
    true = {}
    for s in seqs:
        hits = [(a, b) for f, a, b in jt.get(s, []) if pat.search(f)]
        if hits: true[s] = hits[0]
    ntrue = len(true)
    best = None
    for (job, mo), ins in inst.items():
        if job != r.job: continue
        hit = 0
        for uid, p0, ln in ins:
            if uid in true and p0 >= 0:
                a, b = true[uid]
                if ov(p0 + 1, p0 + ln, a, b) >= 0.3 * (b - a + 1): hit += 1
        if hit >= max(2, 0.5 * len(ins)):
            row = mt[(mt.job == job) & (mt.motif == mo)].iloc[0]
            if best is None or row.bp_obs > best[1]: best = (mo, row.bp_obs, row.bp_exp, row.nhelix_cov, len(ins))
    res.append(dict(job=r.job, group=g, n_used=r.n_used, n_true=ntrue,
                    recovered=best is not None, best_motif=best[0] if best else '',
                    best_bp_obs=best[1] if best else 0, best_bp_exp=best[2] if best else 0,
                    best_hcov=best[3] if best else 0, best_n=best[4] if best else 0))
sens = pd.DataFrame(res); sens.to_csv(OUT + '_sensitivity.tsv', sep='\t', index=False)
d = sens[sens.n_true >= MINTRUE]
print('detectable positive-control jobs', len(d), 'recovered (any motif)', d.recovered.sum())
for lo, hi in [(3, 5), (6, 10), (11, 20), (21, 40), (41, 1000)]:
    s = d[(d.n_true >= lo) & (d.n_true <= hi)]
    if len(s): print(f'n_true {lo}-{hi}: {s.recovered.mean():.2f} ({s.recovered.sum()}/{len(s)})')
