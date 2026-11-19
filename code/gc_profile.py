#!/usr/bin/env python
"""Per-genome profile used to test the CGG reassignment hypothesis.
For each genome:
- tRNA anticodon repertoire (ARAGORN -t -i), with flags for Arg/Gly/Trp tRNAs that could read CGG
  (anticodons CCG, UCG/TCG)
- codon counts over all GTDB/Prodigal CDSs: CGG, CGA, AGG, AGA, GGN, plus total codons
- GC content of the CDS
usage: gc_profile.py LIST.tsv(genome_id<TAB>fna.gz<TAB>faa.gz) OUT.tsv NPROC
"""
import sys, gzip, subprocess, tempfile, os, re, multiprocessing as mp
from collections import Counter
LST, OUT, NP = sys.argv[1], sys.argv[2], int(sys.argv[3])
COMP = str.maketrans('ACGTN', 'TGCAN')

def rfa(fn):
    d = {}; n = None; b = []
    with gzip.open(fn, 'rt') as f:
        for l in f:
            if l[0] == '>':
                if n: d[n] = ''.join(b)
                n = l[1:].rstrip(); b = []
            else: b.append(l.strip())
    if n: d[n] = ''.join(b)
    return d

def work(row):
    gid, fna, faa = row
    try:
        g = {k.split()[0]: v.upper() for k, v in rfa(fna).items()}
        cod = Counter(); gc = at = 0
        for h in rfa(faa):
            p = h.split(' # '); pid = p[0].strip(); s, e, st = int(p[1]), int(p[2]), int(p[3])
            c = pid.rsplit('_', 1)[0]
            if c not in g: continue
            nt = g[c][s - 1:e]
            if st == -1: nt = nt.translate(COMP)[::-1]
            for i in range(0, len(nt) - 2, 3): cod[nt[i:i + 3]] += 1
            gc += nt.count('G') + nt.count('C'); at += nt.count('A') + nt.count('T')
        with tempfile.NamedTemporaryFile('w', suffix='.fna', delete=False) as tf:
            for k, v in g.items(): tf.write(f'>{k}\n{v}\n')
            tmp = tf.name
        out = subprocess.run(['aragorn', '-t', '-i', '-gc1', '-l', '-w', tmp], capture_output=True, text=True).stdout
        os.unlink(tmp)
        trnas = re.findall(r'tRNA-(\w+|\?\S*)\s+\S+\s+\S+\s+\((\w+)\)', out)
        rep = Counter(f'{aa}-{ac.upper().replace("T","U")}' for aa, ac in trnas)
        tot = sum(cod.values())
        cgg_trna = ';'.join(f'{k}:{v}' for k, v in rep.items() if k.endswith('-CCG') or k.endswith('-UCG'))
        return '\t'.join(map(str, [gid, tot, cod['CGG'], cod['CGA'], cod['CGT'] + cod['CGC'], cod['AGA'], cod['AGG'],
                                   sum(cod['GG' + x] for x in 'ACGT'), round(gc / max(1, gc + at), 4), sum(rep.values()),
                                   cgg_trna or 'none', ';'.join(f'{k}:{v}' for k, v in sorted(rep.items()) if k.split('-')[0] in ('Arg', 'Gly', 'Trp'))]))
    except Exception as ex:
        return f'{gid}\tERR\t{ex}'

if __name__ == '__main__':
    rows = [l.rstrip('\n').split('\t') for l in open(LST)]
    with mp.Pool(NP) as pool, open(OUT, 'w') as fo:
        fo.write('genome\tcodons\tCGG\tCGA\tCGY\tAGA\tAGG\tGGN\tgc_cds\tn_trna\tcgg_trnas\targ_gly_trp_trnas\n')
        for r in pool.imap_unordered(work, rows, chunksize=4):
            fo.write(r + '\n')
