# A conserved CGG-triggered +1 frameshift code in the gut bacterial order TANB77

This is a computational discovery made from public genomes (GTDB r226, 143,614 species).

In the human/animal gut family **CAG-508** (order **TANB77**, Clostridia), the arginine codon **CGG** acts as a
conserved **+1 frameshift site**. The C is skipped, and the following GGA/GGU encodes Gly in the +1 frame. The
phenomenon coincides with two changes in the tRNA set: the loss of every CGG-reading tRNA (anticodons CCG/UCG), and
the gain of a single-copy **tRNA-Gly(UCC) variant with G37**. Automated genetic-code inference tools report these sites
as a spurious "CGG→Gly reassignment", and standard gene annotations truncate many TANB77 genes at CGG.

- **Executive summary (PDF):** [`report/executive_summary.pdf`](report/executive_summary.pdf)
- **Full report:** [`report/full_report.md`](report/full_report.md)
- **Project history:** [`plan.md`](plan.md) has the original riboswitch plan and the pivot;
  [`results_interim.md`](results_interim.md) has the interim results, including the retracted "CGG→Gly" interpretation.

## Key results
| Test | Result |
|---|---|
| +1 frame restores homology at CGG (genome frame test) | 74.5% of 1,284 CAG-508 sites; controls 0.7–1.3% |
| Replication in complete genomes | 23/23 single-contig CAG-508 genomes (19 PacBio HiFi, 6 BioProjects) |
| Shared shift sites | 73 sites in ≥2 genera; one site in 13 genera |
| CGG-reading tRNA | 0/23 complete CAG-508 genomes; outgroup RGIG9151 5/5 (no shifting) |
| tRNA-Gly(UCC) G37 variant | 23/23 complete CAG-508 genomes; 0/300 other Clostridia |

**Status:** computational evidence only. Translation-level (proteomic or ribosome-profiling) confirmation is still
needed; see "Testable predictions" in the PDF.

## Repository layout
- **`code/`** — all analysis scripts (Python 3.11: pyhmmer, pandas, numpy, scipy, matplotlib, reportlab, pyrodigal;
  plus ARAGORN, tRNAscan-SE, HMMER, Infernal).
  - **Genetic-code screen:**
    - `gc_screen.py`: marker codon→HMM-column likelihoods.
    - `gc_infer.py`: per-genome and per-clade calls.
    - `gc_genome.py`: whole-proteome Codetta-style re-implementation.
  - **Checks:** `gc_checks.py`, `gc_comparator.py`, `gc_clade.py`, `gc_sites.py`.
  - **Frameshift tests:** `gc_quad.py`, `gc_frames_genome.py`, `gc_genomewide.py`.
  - **tRNA:** `gc_profile.py`, `gc_glyvar.py`, `gc_trna_gly.py`.
  - **Figures and PDF:** `gc_figdata.py`, `gc_example.py`, `make_figs.py`, `build_pdf.py`.
  - **Abandoned riboswitch phase (kept for transparency):** `hmm_pass.py`, `utr_pass.py`, `make_jobs.py`,
    `run_cmfinder.sh`, `parse_motifs.py`, `eval_controls.py`, `rfam_scan.sh`, `orchestrate.sh`,
    `compass_gate.py`, `show_clusters.py`.
- **`report/`** — the PDF, the full report, the figures (`figs/`), and the key result tables (`data/`).

## Reproducing
1. **Inputs (all public):**
   - GTDB r226: `genomic_files_reps`, `bac120/ar53_marker_genes_reps`, metadata, `bac120_r226.tree`.
   - Pfam-A and TIGRFAM 15 HMMs; Rfam 15.
   - Shulgina & Eddy 2021 supplementary files.
   - Melnykov KACI results (Zenodo 10.5281/zenodo.20837797).
   - 23 near-complete CAG-508 genomes from NCBI (accessions in `report/data/cag508_complete.tsv`).
2. **Scripts:** they expect GTDB data under `~/data/` and a taxonomy map at `~/work/taxmap.tsv`, built from GTDB
   metadata. Each script's docstring gives its usage.
3. **Compute:** analyses were run on a 16–32 vCPU / 64–128 GB Linux VM.
   - Avoid running several tRNAscan-SE or full-Pfam jobs in parallel; they are memory hungry.
