import json, sys, urllib.request
R = 'atlas_release_2026-10-r4'
for cid in sys.argv[1:]:
    d = json.load(urllib.request.urlopen(f'https://compass.arcinstitute.org/api/cluster/{cid}?r={R}', timeout=60))
    g = lambda k: d.get(k)
    print(f"== {cid} adj={g('top_adj_arch_name')} fs={g('top_fs_name')} n={g('n_dedup')} len={g('median_length')} "
          f"GC={g('median_gc_pct'):.0f} cov={g('c_observed')}/{g('c_expected')} hel={g('n_helices')}/{g('n_cov_helices')} "
          f"pk={g('has_pseudoknot')} gu={g('gu_frac')} phylum={g('top_phylum')}({g('top_phylum_frac'):.2f}) "
          f"tax95={g('dominant_taxon_95')} meta={g('pct_metagenomic'):.0f}% 3p={g('arch_3p_only_frac')} "
          f"term_any={g('any_terminator')} coor={g('frac_mem_coor')}")
    print('  ', g('consensus')); print('  ', g('dotbracket'))
