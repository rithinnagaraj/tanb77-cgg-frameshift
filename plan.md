# PIVOT (2026-11-19): GTDB-wide screen for undiscovered genetic codes

**Why we pivoted from riboswitches.** The cheap gate against the COMPASS atlas (Evo2-based, 376M intergenic
regions, Oct 2026) looked at all 45 target Pfams:
- CoA genes: every strong structured RNA is a variant of the known pan motif (RF01749; core GUAUCC…GAAC). These
  turned up upstream of panB/panC/panD/panT/coaX in Clostridia, Myxococcota, Acidobacteria and Actinobacteria.
- Biotin genes: the clusters are BirA DNA operators (TTGTnAAC…GTTnACAA) next to a terminator.
- PLP, heme and quinone genes: nothing directly upstream; the hits were r-protein leaders, terminators or
  neighbouring genes.

A novel positive result was therefore unlikely, and the user's bar requires one.

**New question.** Do some of the ~143k GTDB r226 species use genetic codes that have not yet been reported?
- Prior screens were each limited in scope:
  - Shulgina & Eddy 2021: GenBank as of 2020.
  - Fakih et al., NAR 2026: 42k GTDB r214 bacteria with <=1% contamination.
  - Melnykov, PNAS 2026: "thousands" of assemblies.
  - Parks et al., gTranslate 2026: UGA stop codon only.
- About 100k r226 species reps (new species, MAGs above 1% contamination, archaea) have never been screened for
  sense-codon reassignments.

**Method.**
1. Fast screen: codon -> profile-HMM-column likelihoods over the bac120/ar53 marker genes of every species rep,
   pooled by clade.
2. Confirm candidates with full Codetta on whole genomes, tRNA genes (anticodon and identity elements),
   release-factor presence, phylogenetic monophyly, and contig-level consistency.

**Controls.**
- Positives: the known codes must all be recovered (SR1 UGA→Gly + CGA/CGG→Trp; Enterosoma AGG→Met; Anaerococcus
  CGG→Trp; Peptacetobacter CGG→Gln; RAAP-2 ACA→Asp; Mycoplasma/Eggerthellaceae UGA→Trp).
- False positives: standard-code clades, plus shuffled codon-label nulls.

---

# (Superseded) Research plan: Hunting riboswitches for the "blind-spot" cofactors

## Question
Breaker (Biochemistry 2022; TIBS 2022; Curr Biol 2023) highlights near-universal cofactors that have **no known
riboswitch class**: pyridoxal 5'-phosphate (PLP, vitamin B6), biotin, coenzyme A / pantothenate, heme, lipoate,
and quinones (CoQ / menaquinone). Is this a real biological absence, or have searches not reached far enough?
We run a targeted, gene-centric comparative-genomics search for conserved structured RNAs in the 5' regions of
genes for the biosynthesis, salvage and transport of these cofactors. The search covers tens of thousands of
GTDB r226 genomes and is calibrated with positive and negative controls.

Outcomes, all of which give new knowledge:
- (A) One or more strong riboswitch candidates for a blind-spot cofactor. This is the main hoped-for result.
- (B) New structured cis-regulatory RNA motifs (other than riboswitches) at these loci.
- (C) A calibrated negative result: "at sensitivity S, measured on the known cofactor riboswitches, no riboswitch
  class with more than N representatives exists upstream of these genes in the sampled diversity."

## Data
- GTDB r226 representative genomes, sampled as one per genus plus all representatives from under-sampled phyla
  (target: about 20–30k genomes). The Prodigal gene calls come from GTDB.
- Rfam 15 covariance models (to mask known RNAs), and Pfam / TIGRFAM / KOfam HMMs for the target gene families.

## Target gene families (initial)
- PLP: PdxS, PdxT, PdxA, PdxJ, PdxH, PdxK/PdxY, PdxB, PLPBP (YggS), plus B6 transporters (ECF S-units).
- Biotin: BioB, BioF, BioA, BioD, BioW, BioI, BioC, BioH, BioY (ECF S-unit), BioN/BioM.
- CoA / pantothenate: PanB, PanC, PanD, PanE, PanF, PanT, CoaA (PanK I/II/III), CoaBC, CoaD, CoaE.
- Heme: HemA, HemL, HemB, HemC, HemD, HemE, HemN/HemF, HemY/HemG, HemH, ChdC/HemQ, AhbA-D, heme uptake
  (ChuA/HmuTUV), efflux (HrtAB), and degradation (IsdG/HmoB/ChuS).
- Lipoate: LipA, LipB, LplA, LipM, LipL.
- Quinones: UbiA–UbiX, MenA–MenH, MqnA–MqnE.

## Pipeline
1. HMMER search of the proteomes for each gene family. Extract the upstream region of the first gene of each
   predicted operon (−400 to +60 nt around the start codon, bounded by the previous gene).
2. Mask known RNAs with Rfam cmscan, keeping the hits for the control analyses.
3. Cluster the UTR sets by taxon and sequence similarity (MMseqs2 / BLASTN). Run CMfinder 0.4.1.18 on each
   cluster. Expand the motifs iteratively with cmsearch over all intergenic regions, then refine them.
4. Score each motif:
   - R-scape covariation (APC G-test, E < 0.05), the number of significantly covarying pairs, and their power.
   - Consistent 5'-UTR location, gene-family enrichment, a multi-stem architecture, conserved unpaired
     nucleotides, and an expression platform (intrinsic terminator or Shine–Dalgarno sequestration).
5. Run the same pipeline, blind, on control gene families.

## Controls (method validation)
- **Positive controls** (sensitivity): genes known to carry riboswitches, run blind with Rfam labels hidden.
  These are thiC (TPP), ribB/ribD (FMN), btuB/cob (AdoCbl), folT (THF), nadA (NAD+), moaA (Moco), metK (SAM),
  lysC (lysine), gcvT (glycine), crcB (fluoride), queC (preQ1) and mgtE (Mg2+). We will report the recovery
  rate as a function of representative count and phylogenetic spread.
- **Negative controls** (false-positive rate):
  - SISSIz-shuffled versions of every candidate alignment, passed through CMfinder and R-scape.
  - UTRs of housekeeping gene families with no expected riboswitch (e.g., gyrB, recA, rpoB-adjacent
    non-leader families).
  - Synthetic spike-in tests, in which known riboswitch sequences are inserted into random UTR sets at known
    prevalence.
- **False flags to exclude**:
  - Palindromic DNA operators. The BirA/bioO, BioR, PdxR and HrtR/Irr sites are inverted repeats that mimic
    hairpins. Discriminators: multi-stem architecture, G·U-inclusive covariation (meaningful only in RNA), and
    matches to RegPrecise motifs.
  - uORFs and coding regions (RNAcode, BLASTX), REP/BIME and other mobile repeats, ribosomal-protein leaders,
    T-boxes, and lone terminators.

## Verification ladder (increasing confidence)
- L0 Pipeline recovers spiked-in known motifs.
- L1 Blind recovery of the known cofactor riboswitches gives the sensitivity.
- L2 The shuffled and housekeeping nulls give the FDR.
- L3 Candidate covariation: at least 2 stems with significantly covarying pairs, including G·U-inclusive
  covariation, and an E-value from R-scape.
- L4 Held-out genomes: a CM trained on half the genomes predicts occurrences in the other half, and these
  occurrences are again upstream of the same gene family (enrichment p-value against random intergenic regions).
- L5 Phylogenetic breadth and independence: the motif is present in at least 2 orders, and its covariation is
  not driven by near-identical sequences.
- L6 Alternative explanations are excluded: coding, TF operator, mobile element, leader peptide.
- L7 Thermodynamic plausibility: MFE z-score against shuffles, and an alternative-conformation (switch)
  analysis that overlaps the aptamer with the terminator / anti-terminator.
- L8 Novelty: no hit to Rfam 15 (cmscan), none to the published motif supplements (Weinberg 2007/2010/2017,
  Stav 2019, Brewer 2021, Narunsky 2024, Fremin 2021), and a literature check.
- L9 Independent adversarial review of the plan and of the results.

## Testable predictions (if a candidate is found)
- In-line probing or SHAPE shows ligand-dependent structural modulation at specified nucleotides for the
  predicted ligand (e.g., PLP, pyridoxal, or pyridoxamine for pdx-associated motifs), with no response to
  close analogues.
- A reporter fusion in a host carrying the motif responds to the cofactor's availability (OFF or ON, as
  predicted from the expression platform).
- Prevalence predictions for genomes not yet sampled.

## Pivots
- If the blind-spot genes yield nothing, extend the search to Breaker's other blind spots: phosphate
  homeostasis, sugar catabolism, and fatty-acid / lipid genes.
- Otherwise, report outcome (C) together with the motifs from outcome (B).

## Revisions after adversarial plan review (independent subagent)
1. **Check novelty before computing.** Search COMPASS (Arc Institute Evo2 intergenic atlas, Oct 2026), RNAcentral
   and the ZWD supplements. Known motifs at the target loci count as positive controls or prior art, not
   discoveries: pan (RF01749), Zeta-pan (RF03111), hemB/anti-hemB (RF02837), serC (RF00517; PLP did not bind it),
   and the sucA-type motifs.
2. **Use all GTDB r226 species representatives** at the target loci, streaming the genomes and keeping only the
   extracted regions. Group sequences by gene family × clade, not by BLASTN clustering.
3. **Matched null for the FDR:** about 200 random gene families matched on copy number, phylogenetic spread, GC
   content and operon-head status, plus dinucleotide-shuffled and reversed upstream regions.
4. **Add rare classes to the positive controls** (preQ1-III, NiCo, guanidine-II/III, 2'-dG, HMP-PP, NAD+-II, Li+).
   Report a sensitivity curve against representative count and divergence.
5. **Test covariation without the structure-biased alignment.** Re-run R-scape on MAFFT (structure-free) and on
   CaCoFold alignments, and report the power-adjusted number of covarying pairs. Require several independent
   clusters below ~85% identity, each contributing covariation.
6. **Operator discriminators:**
   - Half-site sequence symmetry.
   - Position relative to the predicted promoter (−10/−35).
   - Anti-co-occurrence with the corresponding TF gene (BirA-DBD, BioR, BioQ, PdxR, Irr, HrtR).
7. **Do not discard sub-threshold or variant Rfam hits** at blind-spot genes; flag them for inspection.
8. **L4 holdout by whole clade** (order or phylum), not a random split.
9. **Prioritise CoA, PLP and biotin,** especially in clades that lack the known regulator.
10. **Extra ladder rungs:**
    - Phylogenetic profiling against the pathway (synthesis vs uptake).
    - Conserved nucleotides concentrated in loops and junctions (a pocket signature).
    - CMCompare against all known aptamers.
    - Public Term-seq / dRNA-seq where available.

**Deliverable framing:** the main deliverable is (C) a sensitivity-calibrated census plus (B) the new motifs found.
(A) a riboswitch candidate is upside.

## Compute
Linux VM with about 16–32 vCPU, 64 GB RAM, 400–500 GB SSD, and no GPU. Tools: HMMER, Infernal, CMfinder,
R-scape, MMseqs2, ViennaRNA, BLAST+, RNAcode, TransTermHP, and Python.
