# A conserved CGG-triggered +1 frameshift code in the gut bacterial order TANB77

*Full report accompanying `executive_summary.pdf`. 2026-11-20. Repository: https://github.com/rithinnagaraj/tanb77-cgg-frameshift (code in `code/`).*

## 1. Story of the project, including what failed

1. **Riboswitches (abandoned).** The first plan was a search for riboswitches for "blind-spot" cofactors (PLP, biotin, CoA,
   heme, lipoate, quinones).
   - A cheap novelty gate against the COMPASS atlas (Arc Institute, Oct 2026) showed that every strong structured RNA
     at the target genes was either a known *pan*-motif variant (CoA genes) or a BirA DNA operator (biotin genes).
   - The user requires a novel positive finding, so we pivoted.
2. **Genetic-code screen.** We screened all 143,614 GTDB r226 species representatives with a Codetta-like
   codon→profile-column model over 168 marker HMMs.
   - All known bacterial codes were recovered.
   - One new signal stood out: CGG behaves like Gly in family CAG-508 (TANB77).
3. **First interpretation: "CGG→Gly reassignment".** It survived many checks: codon specificity, contamination, a
   calibrated null, pseudo-replication, and independent KACI calls.
4. **Reinterpretation.** Inspecting alignments showed that the protein downstream of CGG is homologous only in the +1
   frame. CGG is a conserved +1 frameshift site (C-GGA/U), and the "Gly" is encoded by the shifted codon.

## 2. Main results (final numbers)

| Test | Result |
|---|---|
| Genome frame test (marker genes located in genomes; 400 CAG-508 and 200 control genomes) | +1 better at 74.5% of 1,284 CAG-508 CGG sites (median +57.9 bits). CAG-508 GGN controls 1.3% (78,544); CGG in other Clostridia 0.7% (54,478) |
| Position-matched controls (last 10% of annotated gene) | CGG 82% vs GGN 8.6% (edge-effect baseline) |
| Native interior sites (gene GC ≤ 0.333, rel. position 0.5–0.9) | 91.5% (n = 117) |
| Non-shifting native sites | 42 in gene first halves, mostly purB column 3 (an invariant Arg column, 43 genomes): read in frame |
| Insertion comparator: X-GG at Gly columns per GGN codon | CAG-508: CGG 5.0e-4, AGG 3.2e-5, UGG 4.5e-6, i.e. C-specific. AGG/UGG rates are at the other-Clostridia level |
| Next base after CGG, against CAG-508's own GGN usage | A 72% vs 59%; U 18% vs 32% (z≈9.4); C-GGA preferred, consistent with a UCC-anticodon decoder |
| Site sharing | 281 distinct native sites; 73 shared across ≥2 genera; TIGR00580 (Mfd) col 808 in 84 genomes from 13 genera |
| Long-read replication: 23 complete/single-contig CAG-508 genomes from 6 BioProjects (19 PacBio HiFi, 3 Nanopore+Illumina, 1 Illumina) | CGG in gene tails 23–69% (background 5%); +1 continuation median ≈70% vs AGG ≈21% and GGG ≈15% |
| Genome-wide (229 CAG-508 MAGs with ≥90% completeness) | 42.5% of CGGs in the last 15 codons of a gene (other Clostridia 5.3%); 66% with +1 ORF continuation (13%); median 36 such genes per genome, ≈28 above control |
| CGG-reading tRNA (CCG/UCG anticodon) | CAG-508 0.5%, UBA1234 0%, CAG-465 3%, complete CAG-508 0/23. RGIG9151 (TANB77 outgroup) 5/5. Other Clostridia 90% |
| tRNA-Gly(UCC) G37 variant (anticodon loop CUUCCGA vs canonical CUUCCAA) | CAG-508 80% (23/23 complete genomes, single copy), UBA1234 96%, CAG-465 65%, SIG350 2/3. RGIG9151 0/5; other Clostridia 0/300 |
| Phenotype tracks genotype across families | Only families with reader loss and the variant shift; RGIG9151 does not |
| Rare-CGG families GTDB-wide (13) | All keep CGG at Arg sites, except the known Trp clades and CAG-508 |

## 3. Claims, with calibrated confidence

- **C1. Predominant +1 frameshift at CGG in CAG-508.** High.
- **C2. Real, inherited sites.** High: long-read replication plus cross-genus sharing. Read-level SRA mapping was not
  done.
- **C3. Coincides with loss of all CGG-reading tRNAs.** Moderate–high: a single evolutionary contrast, with complete
  genomes concordant.
- **C4. A co-occurring tRNA-Gly(UCC) G37 variant plausibly decodes C-GGA.** Correlative; the mechanism is a
  hypothesis. Within MAGs, per-genome presence is not predictive because tRNA recovery in MAGs is incomplete. In
  complete genomes the variant is present in 23 of 23.
- **C5. Annotations truncate genes, and code-inference tools report a frameshift shadow as "CGG→Gly".** High.
- **C6. A minority of CGGs are read in frame** (for example purB). Moderate.

## 4. Novelty

- **Literature and data checks** are as in the PDF, §3.
- **Prior art that narrows the claim:**
  - Isolated +1 shifting at rare or uncharged codons: E. coli CGG-C/AGG-AGG (Spanjaard & van Duin 1988; FSscan).
  - Recoded E. coli resolving unassigned codons by frameshifting (Ma 2018).
  - Yeast near-cognate codons (Sundararajan 1999).
  - Eukaryotic frameshift codes (Euplotes; Perkinsus mitochondria).
  - The KACI table contains undescribed "CGG=G" calls for CAG-508.
- **What is new:** a clade-wide, conserved, tRNA-linked frameshift code in bacteria. It co-occurs with a dedicated
  tRNA-Gly variant, and it explains a class of spurious code-inference calls.

## 5. Open items and next steps

1. Read-level support: map SRA reads at about 30 shift sites in short-read and HiFi data, and check C-allele
   frequency and strand balance.
2. A frameshift-aware census of all proteins (DIAMOND `-F 15`) in CAG-508 versus matched Clostridia: all shift
   contexts, and the overall indel rate.
3. Selection: whether shift sites are depleted in highly expressed genes, and dN/dS in the +1 frame across genera.
4. tRNA identity elements of the G37 variant (tRNAscan-SE isotype and structure) versus its A37 paralog; the
   co-phylogeny of the variant.
5. Re-test other Clostridial CGG "reassignments" (Gln/Trp/Ala) with the frame test.
6. Experimental: reporter assays with the variant tRNA, and proteomics of enriched TANB77.

## 6. Files

- **Figures:** `figs/fig0–5`.
- **Results tables (VM `~/gc/res/`):**
  - `frames_genome.tsv`
  - `frames_cag508_gc.tsv`
  - `genomewide_annot.tsv`
  - `genomewide_lr.tsv`
  - `glyvar_annot.tsv`
  - `glyvar_lr.tsv`
  - `cag508_cgg_sites_gc.tsv`
  - `null_rare_codons.tsv`
  - `checks.txt`
  - `profile_annot.tsv`
- **Interim notes:** `../results_interim.md`.
