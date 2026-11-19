#!/usr/bin/env python
"""Run tRNAscan-SE on a genome sample and extract the anticodon arms (stem + loop) of the Gly, Arg and Trp
tRNAs, to look for clade-specific tRNA variants that could decode CGG.
usage: gc_trna_gly.py LIST.tsv(genome<TAB>fna.gz<TAB>faa) OUT.tsv NPROC
"""
import sys, gzip, subprocess, tempfile, os, re, multiprocessing as mp
LST, OUT, NP = sys.argv[1], sys.argv[2], int(sys.argv[3])

def work(row):
    gid, fna = row[0], row[1]
    td = tempfile.mkdtemp()
    fa = os.path.join(td, 'g.fna'); ss = os.path.join(td, 'g.ss')
    with gzip.open(fna, 'rt') as f, open(fa, 'w') as o: o.write(f.read())
    subprocess.run(['tRNAscan-SE', '-B', '-q', '--brief', '-f', ss, '-o', os.devnull, fa], capture_output=True)
    out = []
    if os.path.exists(ss):
        txt = open(ss).read().split('\n\n')
        for blk in txt:
            m = re.search(r'Type: (\w+)\s+Anticodon: (\w+) at (\d+)-(\d+)', blk)
            sq = re.search(r'Seq: (\S+)', blk); st = re.search(r'Str: (\S+)', blk)
            sc = re.search(r'Score: ([\d.]+)', blk)
            if not (m and sq and st): continue
            iso, ac, a0 = m.group(1), m.group(2), int(m.group(3))
            if iso not in ('Gly', 'Arg', 'Trp', 'Undet'): continue
            s = sq.group(1).upper(); a = a0 - 1
            loop = s[a - 2:a + 5] if a >= 2 else ''
            arm = s[a - 7:a + 10] if a >= 7 else ''
            stru = st.group(1)[a - 7:a + 10] if a >= 7 else ''
            out.append('\t'.join([gid, iso, ac, sc.group(1) if sc else '', loop, arm, stru, s]))
    subprocess.run(['rm', '-rf', td])
    return out

if __name__ == '__main__':
    rows = [l.rstrip('\n').split('\t') for l in open(LST)]
    with mp.Pool(NP) as pool, open(OUT, 'w') as fo:
        fo.write('genome\tisotype\tanticodon\tscore\tac_loop\tac_arm\tac_arm_str\tseq\n')
        for res in pool.imap_unordered(work, rows):
            fo.write(''.join(r + '\n' for r in res))
