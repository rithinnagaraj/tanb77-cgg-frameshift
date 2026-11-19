#!/usr/bin/env python
"""Cheap novelty / feasibility gate against the COMPASS atlas (public website API, highlighted scope).

For each target Pfam the script lists the COMPASS clusters whose loci are strongly associated with that
gene (g0_pct, the % of loci) and that show covariation with no Rfam annotation. For each such cluster it
fetches the consensus structure and alignment and flags signs of a DNA operator: a long, near-perfect
inverted repeat in the consensus (palindrome score) and few helices.
"""
import json, time, sys, urllib.request, re, collections

B = 'https://compass.arcinstitute.org/'
R = 'atlas_release_2026-10-r4'
PF = {
 'PLP': {'PF01680': 'PdxS', 'PF01174': 'PdxT', 'PF04166': 'PdxA', 'PF03740': 'PdxJ', 'PF01243': 'PdxH',
         'PF10590': 'PdxH_C', 'PF08543': 'PdxK'},
 'BIOTIN': {'PF06968': 'BioB_BATS', 'PF13500': 'BioD', 'PF03744': 'BioW', 'PF02632': 'BioY', 'PF03099': 'BirA/LplA'},
 'COA': {'PF02548': 'PanB', 'PF02569': 'PanC', 'PF02261': 'PanD', 'PF02558': 'PanE', 'PF08546': 'PanE_C',
         'PF04127': 'CoaB(DFP)', 'PF02441': 'CoaC/UbiX', 'PF03630': 'CoaW', 'PF03309': 'CoaX', 'PF01121': 'CoaE'},
 'HEME': {'PF05201': 'HemA', 'PF00490': 'HemB', 'PF01379': 'HemC', 'PF02602': 'HemD', 'PF01208': 'HemE',
          'PF06969': 'HemN', 'PF01218': 'HemF', 'PF12724': 'HemG', 'PF00762': 'HemH', 'PF06778': 'ChdC'},
 'QUINONE': {'PF01977': 'UbiD', 'PF02621': 'MqnA/D', 'PF01040': 'UbiA/MenA'},
}

def post(path, body):
    req = urllib.request.Request(B + path, data=json.dumps(body).encode(), headers={'Content-Type': 'application/json'})
    return json.load(urllib.request.urlopen(req, timeout=60))

def get(path):
    return json.load(urllib.request.urlopen(B + path, timeout=60))

def rc(s): return s[::-1].translate(str.maketrans('ACGU', 'UGCA'))

def palindrome_score(seq):
    """Longest near-perfect inverted repeat (arm length, <=1 mismatch, loop 0-20 nt) in the consensus."""
    s = seq.upper().replace('T', 'U'); best = 0
    for L in range(6, 21):
        for i in range(0, len(s) - 2 * L):
            arm = s[i:i + L]
            for loop in range(0, 21):
                j = i + L + loop
                if j + L > len(s): break
                mm = sum(1 for a, b in zip(arm, rc(s[j:j + L])) if a != b)
                if mm <= 1: best = max(best, L); break
    return best

out = []
for grp, pfs in PF.items():
    for pf, name in pfs.items():
        q = {'q': '', 'mode': 'all', 'conds': [{'t': 'gene', 'kind': 'pfam', 'terms': [pf], 'min': 0}]}
        try: res = post('api/search', {'query': q, 'sort': 'c_observed', 'desc': True, 'cols': [], 'limit': 2000})
        except Exception as ex:
            print('ERR', pf, ex, file=sys.stderr); continue
        rows = res['rows']
        for r in rows:
            g = r.get('g0_pct') or 0; c = r.get('c_observed') or 0
            if g < 40 or c < 8 or r.get('top_rfam'): continue
            cid = r['cluster_id']
            try: st = get(f'api/structure/{cid}?r={R}')
            except Exception: st = {}
            cons = st.get('consensus', ''); db = st.get('dotbracket', '')
            out.append(dict(group=grp, pfam=pf, gene=name, cluster=cid, pct_loci=g, cov=c, nseq=r.get('n_dedup'),
                            len=r.get('median_length'), gc=r.get('median_gc_pct'), term=r.get('frac_mem_coor'),
                            adj=r.get('top_adj_arch_name'), tier=r.get('tier_label'), helices=st.get('helices'),
                            pairs=st.get('pairs'), stcov=st.get('covarying_pairs'), pal=palindrome_score(cons) if cons else None,
                            dotbracket=db, consensus=cons))
            time.sleep(0.2)
        print(grp, pf, name, 'clusters', res.get('total'), 'candidates so far', len(out), file=sys.stderr)
        time.sleep(0.5)
import csv
with open(sys.argv[1], 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(out[0].keys()), delimiter='\t'); w.writeheader(); w.writerows(out)
