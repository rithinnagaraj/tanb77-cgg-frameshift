#!/usr/bin/env python
"""Whole-proteome, Codetta-style genetic code inference for one genome (an independent re-implementation).

Inputs: the genome FASTA (.fna.gz), the GTDB Prodigal proteins (.faa.gz, whose headers carry CDS coordinates),
and a Pfam HMM file.

Method:
1. hmmsearch all Pfam profiles against the proteome (domain i-Evalue <= 1e-10). Profiles on Codetta's
   exclusion list (transposons, viral, etc.) are skipped.
2. Map every aligned match-state column to the codon in the CDS that encodes the residue.
3. For each codon c and amino acid a, logL[c,a] = sum over columns of log e_col(a). Columns are capped
   Codetta-style: each (profile, column) counts at most once per codon.
4. The decoding probability is P(a | c), proportional to exp(logL[c,a]) (uniform prior). A codon is called
   when P >= 0.9999.
Optional: --exclude-markers FILE lists marker Pfam/TIGR accessions to drop, so that the evidence is
independent of the marker screen.

usage: gc_genome.py GENOME.fna.gz PROTEINS.faa.gz PFAM.hmm OUT_PREFIX [--exclude-markers FILE] [--bad FILE] [--cpus N]
"""
import sys, gzip, argparse, numpy as np, pyhmmer, json

ap = argparse.ArgumentParser()
ap.add_argument('fna'); ap.add_argument('faa'); ap.add_argument('hmm'); ap.add_argument('out')
ap.add_argument('--exclude-markers'); ap.add_argument('--bad'); ap.add_argument('--cpus', type=int, default=4)
ap.add_argument('--evalue', type=float, default=1e-10)
ap.add_argument('--dump-codons', default='CGG', help='comma-separated codons whose aligned sites are written to OUT.sites.tsv')
A = ap.parse_args()
AA = 'ACDEFGHIKLMNPQRSTVWY'
B = 'TCAG'; COD = [a + b + c for a in B for b in B for c in B]; CIDX = {c: i for i, c in enumerate(COD)}
STD = dict(zip(COD, 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
COMP = str.maketrans('ACGTN', 'TGCAN')

def rfa(fn):
    op = gzip.open if fn.endswith('.gz') else open
    d = {}; n = None; buf = []
    with op(fn, 'rt') as f:
        for l in f:
            if l[0] == '>':
                if n: d[n] = ''.join(buf)
                n = l[1:].rstrip('\n'); buf = []
            else: buf.append(l.strip())
    if n: d[n] = ''.join(buf)
    return d

genome = {k.split()[0]: v.upper() for k, v in rfa(A.fna).items()}
prots = rfa(A.faa)
cds = {}; pseq = []
alpha = pyhmmer.easel.Alphabet.amino()
for h, s in prots.items():
    p = h.split(' # '); pid = p[0].strip(); st, en, strand = int(p[1]), int(p[2]), int(p[3])
    contig = pid.rsplit('_', 1)[0]
    if contig not in genome: continue
    nt = genome[contig][st - 1:en]
    if strand == -1: nt = nt.translate(COMP)[::-1]
    cds[pid] = nt
    pseq.append(pyhmmer.easel.TextSequence(name=pid.encode(), sequence=s.replace('*', '').replace('X', 'A')).digitize(alpha))
blk = pyhmmer.easel.DigitalSequenceBlock(alpha, pseq)

skip = set()
for f in (A.exclude_markers, A.bad):
    if f: skip |= set(x.split('.')[0] for x in open(f).read().split())

hmms = []
with pyhmmer.plan7.HMMFile(A.hmm) as f:
    for h in f:
        acc = h.accession or h.name; acc = acc.decode() if isinstance(acc, bytes) else acc
        if acc.split('.')[0] in skip: continue
        hmms.append(h)
lem = {}
logL = np.zeros((64, 20)); n = np.zeros(64, int); seen = [set() for _ in range(64)]; Z = np.zeros(20)
DUMP = set(A.dump_codons.split(',')); sites = []
nhits = 0
for hits in pyhmmer.hmmsearch(hmms, blk, cpus=A.cpus, E=1e-5, domE=A.evalue):
    hm = hits.query
    qn = hm.accession or hm.name; qn = qn.decode() if isinstance(qn, bytes) else qn
    if qn not in lem:
        em = np.clip(np.array(hm.match_emissions)[1:, :20], 1e-6, 1); lem[qn] = np.log(em)
    L = lem[qn]
    for hit in hits:
        tname = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
        nt = cds.get(tname)
        if nt is None: continue
        for dom in hit.domains:
            if dom.i_evalue > A.evalue: continue
            ali = dom.alignment; nhits += 1
            k = ali.hmm_from - 1; r = ali.target_from - 1
            for hc, tc in zip(ali.hmm_sequence, ali.target_sequence):
                if hc == '.':                  # insert state in target
                    r += 1; continue
                if tc == '-':                  # deletion in target
                    k += 1; continue
                c = CIDX.get(nt[3 * r: 3 * r + 3])
                if c is not None and COD[c] in DUMP:
                    e = np.exp(L[k])
                    sites.append(f'{tname}\t{tname.rsplit("_", 1)[0]}\t{r + 1}\t{COD[c]}\t{qn}\t{k + 1}\t{AA[int(e.argmax())]}\t{e[AA.index("G")]:.3f}\t{e[AA.index("R")]:.3f}\t{tc.upper()}')
                if c is not None and (qn, k) not in seen[c]:
                    seen[c].add((qn, k)); logL[c] += L[k]; n[c] += 1; Z += np.exp(L[k])
                k += 1; r += 1
res = {}
# Codetta likelihood: P(column k | codon decodes a) = e_k(a) / sum_j e_j(a), over all observed columns j;
# the 21st "nonspecific" model draws columns uniformly.
logZ = np.log(Z); ntot = n.sum()
for c in range(64):
    m = n[c]
    if m == 0:
        s = np.full(21, np.log(1 / 21))
    else:
        s = np.empty(21)
        s[:20] = np.log(1 / 21) + logL[c] - m * logZ
        s[20] = np.log(1 / 21) - m * np.log(ntot)
    p = np.exp(s - s.max()); p /= p.sum(); a = int(p.argmax())
    lab = (AA + '?')[a]
    res[COD[c]] = dict(std=STD[COD[c]], call=lab if (p[a] >= 0.9999 and a < 20) else '?', p=float(p[a]), n=int(m),
                       top3=[((AA + '*')[i], float(p[i])) for i in np.argsort(-p)[:3]])
code = ''.join(res[c]['call'] for c in COD)
json.dump(dict(code=code, domains=nhits, codons=res), open(A.out + '.json', 'w'), indent=1)
with open(A.out + '.sites.tsv', 'w') as f:
    f.write('protein\tcontig\taa_pos\tcodon\tprofile\tcol\tcol_consensus\te_G\te_R\ttranslated\n')
    f.write(''.join(s + '\n' for s in sites))
# per-contig composition (for contamination analysis): length, GC, CGG count and Arg-codon usage over CDS
from collections import Counter
cc = {}
for pid, ntc in cds.items():
    ctg = pid.rsplit('_', 1)[0]; d = cc.setdefault(ctg, Counter())
    for i in range(0, len(ntc) - 2, 3): d[ntc[i:i + 3]] += 1
with open(A.out + '.contigs.tsv', 'w') as f:
    f.write('contig\tlength\tgc\tcodons\tCGG\tCGN\tAGR\n')
    for ctg, d in cc.items():
        s = genome[ctg]; gcf = (s.count('G') + s.count('C')) / max(1, len(s))
        f.write(f'{ctg}\t{len(s)}\t{gcf:.4f}\t{sum(d.values())}\t{d["CGG"]}\t{sum(d["CG" + x] for x in "ACGT")}\t{d["AGA"] + d["AGG"]}\n')
print('domains', nhits)
print('Genetic code:', code)
print('Standard    :', ''.join(STD[c] for c in COD))
for c in COD:
    r = res[c]
    if r['call'] not in ('?', r['std']): print('NONSTANDARD', c, r['std'], '->', r['call'], 'n=', r['n'], 'p=', round(r['p'], 6))
