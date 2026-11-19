#!/usr/bin/env python
"""Collect CMfinder motifs and their R-scape statistics into one table.

For each motif the table records:
- nseq, alignment length, average identity, number of base pairs, and number of helices (>= 3 bp)
- the number of covarying base pairs observed vs. expected (R-scape power), and the number of helices
  carrying >= 1 significant covarying pair
- the fraction of the job's sequences that contain the motif
- where the motif sits relative to the start codon (median distance of its 3' end upstream of the ATG)
- the number of distinct genomes and classes, and the max/median pairwise identity of the motif sequences
It also writes motifs.fa, holding every motif instance, for the later Rfam / novelty scan.

usage: parse_motifs.py JOBDIR UTR_TSV OUT_TSV OUT_FA
"""
import sys, os, re, glob, statistics
import pandas as pd

JD, UTSV, OUT, OFA = sys.argv[1:5]
u = pd.read_csv(UTSV, sep='\t', usecols=['uid', 'up_len', 'genome']).set_index('uid')
UP, IN = 400, 20

def read_sto(fn):
    seqs = {}; ss = ''
    for l in open(fn):
        if l.startswith('#=GC SS_cons'): ss += l.split()[-1]
        elif l.startswith('#') or l.startswith('//') or not l.strip(): continue
        else:
            n, s = l.split()[:2]; seqs[n] = seqs.get(n, '') + s
    return seqs, ss

def pairs_of(ss):
    st = {}; pr = []
    br = {'<': '>', '(': ')', '[': ']', '{': '}'}
    rev = {v: k for k, v in br.items()}
    for i, c in enumerate(ss):
        if c in br: st.setdefault(c, []).append(i)
        elif c in rev and st.get(rev[c]): pr.append((st[rev[c]].pop(), i))
    return sorted(pr)

def helices(pr):
    hs = []; cur = []
    for p in pr:
        if cur and p[0] == cur[-1][0] + 1 and p[1] == cur[-1][1] - 1: cur.append(p)
        else:
            if cur: hs.append(cur)
            cur = [p]
    if cur: hs.append(cur)
    # merge helices separated by small bulges (<=2 nt either side)
    hs.sort()
    merged = []
    for h in hs:
        if merged and 0 < h[0][0] - merged[-1][-1][0] <= 3 and 0 < merged[-1][-1][1] - h[0][1] <= 3:
            merged[-1].extend(h)
        else: merged.append(h)
    return merged

def parse_rscape(fn):
    d = {'bp_obs': 0, 'bp_exp': 0.0, 'bp_exp_sd': 0.0, 'nbp_rs': 0, 'cov_pairs': [], 'avgid': None}
    if not os.path.exists(fn): return d
    inpow = False
    for l in open(fn):
        m = re.search(r'avgid ([\d.]+)', l)
        if m and l.startswith('# MSA'): d['avgid'] = float(m.group(1))
        if l.startswith('# BPAIRS observed to covary'): d['bp_obs'] = int(l.split()[-1])
        elif l.startswith('# BPAIRS expected to covary'):
            p = l.split(); d['bp_exp'] = float(p[5]); d['bp_exp_sd'] = float(p[7])
        elif l.startswith('# BPAIRS '):
            try: d['nbp_rs'] = int(l.split()[2])
            except: pass
        if l.startswith('# covary'): inpow = True; continue
        if inpow:
            if l.startswith('#-'): continue
            if l.startswith('#') or not l.strip(): inpow = False; continue
            p = l.split()
            if len(p) >= 5 and p[0] == '*':
                d['cov_pairs'].append((int(p[1]), int(p[2]), int(p[3])))
    return d

def ident(a, b):
    m = n = 0
    for x, y in zip(a, b):
        if x in '.-' and y in '.-': continue
        n += 1; m += (x == y)
    return m / n if n else 0

rows = []
fa = open(OFA, 'w')
for jdir in sorted(glob.glob(os.path.join(JD, '*.cmf'))):
    job = os.path.basename(jdir)[:-4]
    jobseq = {}; n = None
    for l in open(os.path.join(jdir, job + '.fa')):
        if l[0] == '>': n = l[1:].strip(); jobseq[n] = ''
        else: jobseq[n] += l.strip()
    lst = os.path.join(jdir, 'motifs.list')
    if not os.path.exists(lst): continue
    for m in open(lst).read().split():
        fn = os.path.join(jdir, m)
        if not os.path.exists(fn): continue
        seqs, ss = read_sto(fn)
        if not seqs: continue
        pr = pairs_of(ss); hs = [h for h in helices(pr) if len(h) >= 3]
        rs = parse_rscape(os.path.join(jdir, 'rs_' + m, 'out.txt'))
        covset = set((a - 1, b - 1) for a, b, _ in rs['cov_pairs'])
        # R-scape coordinates refer to alignment columns (1-based), matching SS_cons columns
        hcov = sum(1 for h in hs if any(p in covset for p in h))
        dists = []; genomes = set()
        for name, s in seqs.items():
            uid = name.split('/')[0]
            raw = s.replace('-', '').replace('.', '').upper().replace('U', 'T')
            js = jobseq.get(uid, '')
            pos = js.find(raw)
            if uid in u.index:
                genomes.add(u.at[uid, 'genome'])
                upl = min(int(u.at[uid, 'up_len']), UP)
                if pos >= 0: dists.append(upl - (pos + len(raw)))  # nt between motif 3' end and ATG (neg = overlaps CDS)
            fa.write(f'>{job}|{m.split(".fa.")[-1]}|{uid}|{pos}|{len(raw)}\n{raw}\n')
        names = list(seqs)
        ids = [ident(seqs[names[i]], seqs[names[j]]) for i in range(len(names)) for j in range(i + 1, min(len(names), i + 15))]
        rows.append(dict(job=job, ko=job.split('__')[0], clade=job.split('__', 1)[1], motif=m.split('.fa.')[-1],
                         nseq=len(seqs), njob=len(jobseq), frac=len(seqs) / max(1, len(jobseq)),
                         alen=len(ss), nbp=len(pr), nhelix=len(hs), nhelix_cov=hcov,
                         bp_obs=rs['bp_obs'], bp_exp=rs['bp_exp'], bp_exp_sd=rs['bp_exp_sd'],
                         avgid=rs['avgid'], med_pid=statistics.median(ids) if ids else None,
                         dist_atg_med=statistics.median(dists) if dists else None,
                         ngenomes=len(genomes)))
pd.DataFrame(rows).to_csv(OUT, sep='\t', index=False)
print('motifs', len(rows))
