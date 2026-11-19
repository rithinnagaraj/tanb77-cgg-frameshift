#!/usr/bin/env python
"""Fast survey of the tRNA-Gly(UCC) anticodon-loop variant with ARAGORN (-fo, FASTA output).
For each genome, extract every tRNA-Gly(tcc) sequence. Classify the anticodon loop by locating 'TTCC' + 2 nt:
canonical 'CTTCCAA' (A37) versus variant 'CTTCCGA' (G37). Also record whether any CCG/TCG-anticodon tRNA exists.
usage: gc_glyvar.py LIST.tsv OUT NPROC
"""
import sys, gzip, subprocess, tempfile, os, re, multiprocessing as mp
def work(row):
    gid, fna = row[0], row[1]
    with tempfile.NamedTemporaryFile('w', suffix='.fna', delete=False) as tf:
        with gzip.open(fna, 'rt') as f: tf.write(f.read())
        tmp = tf.name
    out = subprocess.run(['aragorn', '-t', '-gc1', '-l', '-fo', tmp], capture_output=True, text=True).stdout
    os.unlink(tmp)
    recs = re.findall(r'>(\S+)[^\n]*\n([acgtn\n]+)', out)
    gly = [s.replace('\n', '').upper() for h, s in recs if 'tRNA-Gly(tcc)' in h]
    var = sum(1 for s in gly if 'CTTCCGA' in s); can = sum(1 for s in gly if 'CTTCCAA' in s)
    cgg = sum(1 for h, s in recs if re.search(r'tRNA-\w+\((ccg|tcg)\)', h))
    return f'{gid}\t{len(gly)}\t{var}\t{can}\t{cgg}'
if __name__ == '__main__':
    rows = [l.rstrip('\n').split('\t') for l in open(sys.argv[1])]
    with mp.Pool(int(sys.argv[3])) as p, open(sys.argv[2], 'w') as fo:
        fo.write('genome\tn_glyUCC\tvariant_G37\tcanonical_A37\tn_CGG_reader\n')
        for r in p.imap_unordered(work, rows): fo.write(r + '\n')
