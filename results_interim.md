# Interim results: a CGG→glycine codon reassignment in the gut bacterial family CAG-508 (order TANB77, Clostridia)

## Context
- **Goal:** a novel, publishable discovery, verified computationally.
- **Pivot:** the first direction (riboswitches for cofactors) was abandoned after a COMPASS-atlas gate, which found only
  pan-motif variants and BirA operators. We pivoted to a GTDB r226-wide screen for undiscovered genetic codes.

## Method (screen)
- **Data:** the 120 bac120 + 53 ar53 GTDB marker genes (168 distinct HMMs) for all 143,614 GTDB r226 species
  representatives.
- **Alignment:** each marker protein was aligned to its Pfam/TIGRFAM profile with pyhmmer, and each match column was
  mapped back to its codon in the marker CDS.
- **Per (genome, codon):** the sum of log emission probabilities over *consensus* columns (max emission ≥ 0.5),
  with Codetta's normalisation:

      P(col | aa) = e(aa) / Σ_cols e(aa)

  plus a nonspecific model. The posterior is computed over 20 aa + nonspecific; a codon is called at P ≥ 0.9999
  with N ≥ 5 columns.
- **Pooling:** evidence was pooled per genus and per family.
- **Independent confirmation:** our own whole-proteome re-implementation of Codetta (all Pfam-A domains, i-E ≤ 1e-10,
  Codetta's bad-Pfam list, each column counted once per codon, bac120 marker Pfams excluded so the evidence is
  independent of the screen).

## Controls recovered (positive controls)
- Genome-level non-standard calls (143k genomes): UGA→Trp 427, UGA→Gly 204, CGA→Trp 56, CGG→Trp 4, AGG→Met 1, and
  one ATG→Leu outlier.
- These correspond to the known codes: Mycoplasmatales and others for code 4; BD1-5/Absconditabacterales for
  code 25/37; Absconditabacterales CGA/CGG→Trp; Enterosoma AGG→Met.
- Pooled genus/family calls also recover RFN20/Onthovivens CGG→Trp and Absconditabacterales CGA→Trp (with many
  individually supporting genomes).
- Anaerococcus (CGG→Trp, a known code) shows Trp as the top CGG posterior in our whole-genome implementation, but
  with only 20 columns (p = 0.49, i.e. underpowered).

## Artefacts identified and handled
- Pooled-only calls with **zero** individually supporting genomes appeared for ATG→C/L/Y/I and AGY→G/T/C, mostly in
  archaea and high-GC clades.
- These are model biases from approximating the Codetta normaliser with marker-HMM-average emissions.
- Filter: require individual-genome support or independent confirmation.

## The candidate: CGG (Arg) → Gly in family CAG-508 (830 species, order TANB77, Clostridia; human/animal gut)

1. **Pooled marker evidence.** CGG is called G with posterior ≈ 1 (391 consensus columns). CGG aligned to
   Gly-consensus columns 332 times and Arg-consensus columns 47 times. Five TANB77 genera are called G individually
   at genus level (CAG-269 with 322 species, JAAWPK01, RGIG8482, UMGS1663, UMGS1781). No other clade in the whole
   GTDB screen gave a CGG→Gly call.

2. **Codon specificity (rules out sequencing/assembly errors).** Fraction of each codon aligned to Gly-consensus
   columns in CAG-508:

   | Codon | Fraction at Gly columns | Columns |
   |---|---|---|
   | CGG | 0.849 | 391 |
   | CGA | 0.0011 | 2,850 |
   | CGU | 0.0006 | 21,663 |
   | CGC | 0.0016 | 1,254 |
   | AGA | 0.0018 | 127,950 |
   | AGG | 0.0023 | 9,213 |

   For CGG in other Clostridia orders the fraction is 0.0009–0.0101 (36k–157k columns). G→C errors at codon
   position 1 would also put CGA/CGU at Gly sites, and they do not.

3. **Site conservation.** Of 1,250 CGG occurrences in CAG-508 marker genes (all match columns), the reference
   Clostridia residue at the same column is Gly-majority at 800, Arg-majority at 91, and other at 359. 474 sites
   are ≥95% Gly in the reference; 61 are ≥95% Arg.

4. **Contamination explains the Arg signal.** CGG at Arg-conserved sites (n = 63) is on marker gene copies with
   median GC 0.401, and 98.4% of them are above the 95th percentile of CAG-508 marker GC (0.333). CGG at
   Gly-conserved sites (n = 566) is on gene copies with median GC 0.304 (background median 0.290); 15% are above
   the 95th percentile. Restricted to native-GC copies (979 sites), the reference-majority residue is Gly at 71%
   and Arg at under 2%, and every CAG-508 genus is Gly-dominant.

5. **Prior art in the 2021 data.** In the Shulgina & Eddy 2021 Codetta results (supp1), 68 CAG-508 genomes have
   CGG called '?' in 58 and 'R' in 10. The 10 'R' calls are scattered across Gly-dominant genera, which is
   consistent with MAG contamination: native CGG usage is about 1.7/10k marker codons versus 57–99 in other
   Clostridia, so a 1–2% contaminant dominates the vote. Their supp2 (all candidate codes) has **no** CGG→Gly
   entry.

6. **Whole-genome independent inference** (marker Pfams excluded in run 2). In 5 TANB77 genus representatives,
   the CGG top call is G in all 5 (2 at P ≥ 0.9999; others 0.97–0.9996). Every other codon is standard. A larger
   run of 90 genomes (48 CAG-508, 23 sister TANB77, 12 other Clostridia, and the 10 Codetta-'R' genomes) was in
   progress when the VM became unreachable: 50/90 finished; early CAG-508 results mostly G-leaning, with few
   columns per genome.

7. **tRNA repertoire** (ARAGORN, 1,228 genomes; tRNAscan-SE on 5).
   - No tRNA with a CCG or UCG anticodon (the only tRNAs that can read CGG by standard pairing) in 0/83 high-quality
     CAG-508 genomes (2/830 overall, likely contamination). Other Clostridia: CCG tRNA present in 73%.
   - Sister TANB77 families have also mostly lost it (7%); CGG is rare across the whole low-GC (~0.28) order.
   - CAG-508 Arg tRNAs: ACG (→ICG, reading CGU/C/A), CCU (AGG), UCU (AGA). Gly tRNAs: GCC and UCC (×2–3),
     including an unusual, conserved tRNA-Gly(UCC) variant with G37 and an anticodon-stem mismatch.
   - The decoding tRNA for CGG→Gly is unknown. Near-cognate candidates: tRNA-Gly(UCC) with a codon-position-1
     C·C mismatch, or an unannotated tRNA.

8. **Novelty check.**
   - None of these report CGG→Gly: the NCBI tables; Shulgina & Eddy 2021 (main + supp2); the Nov 2026 bioRxiv list
     of bacterial codes 4/25/34–37; Fakih 2026 (ACA→Asp); Melnykov PNAS 2026 (Asgard sense reassignment, sense→stop);
     Parks 2026 (Eggerthellaceae UGA→Trp); Kivenson 2025 (Pyl code); gTranslate.
   - Arg→Gly reassignment is otherwise known only in tunicate mitochondria, at AGA/AGG.
   - Still to check: the contents of Melnykov's Zenodo candidate list.

## UPDATE (after review, Nov 20): reinterpretation — CGG is a +1 frameshift / quadruplet (CGGR) site

**Novelty against KACI.** Melnykov's deposited KACI results (Zenodo 20837797, 3.1M assemblies) contain automated
CGG=G calls for 863 CAG-508 assemblies. 98.6% of all GTDB-mapped KACI CGG→G calls fall in CAG-508; the rest are in
TANB77 sister families. His paper does **not** describe them; it explicitly calls Clostridial non-standard CGG decoding
"a subject of a separate investigation". So this is independent-method corroboration, and the finding itself is
unreported.

**Review checks on the marker screen, all passed:**
- Rate ratio of CGG@Gly per GGG versus CGH@Gly per GGH: 439x (p≈0).
- G versus A: +272 nats.
- One genome per genus (44 genera): P(G) = 1. Leaving out CAG-269: P(G) = 1.
- 281 distinct sites: 157 Gly-majority, 10 Arg-majority.
- Empirical null over 812 rare family×codon combinations: p = 0.0012, and the single higher value is an Ile/Leu
  similarity artefact.
- In 13 families with rare CGG, all others keep CGG at Arg sites (97–100%), including SAR11, Clostridiaceae and
  Tissierellaceae.
- Tree: RGIG9151, the family branching outside the others in TANB77, keeps tRNA-CCG and Arg. CAG-465, SIG350, UBA1234
  and CAG-508 lost tRNA-CCG, consistent with a single loss. The Gly signal is strongest in CAG-508.

**Alignment inspection changes the interpretation.** At the most-hit conserved-Gly sites (TIGR00580 col 808,
TIGR02386 col 696, TIGR00445 col 234), the CAG-508 sequence carries an extra C immediately before a GGA Gly codon.
- Annotated frame: ...CGG A.. — Prodigal reads CGG, then junk, then an early stop.
- The +1 frame (GGA ...) restores exact homology downstream (e.g. GNMNQLRQIAGMR, GEIQHGHMEQVGY, GGVISALALYLKM).
- The same site is shifted in genomes from several genera.

**Hypothesis.** In CAG-508, CGG lacks any cognate tRNA, and CGGR acts as a +1 frameshift / quadruplet site:
- One Gly is inserted, then decoding continues in +1, either through hungry-codon P-site slippage or through
  quadruplet pairing by the conserved G37 variant of tRNA-Gly(UCC): 5'-U34 C35 C36 G37 pairing with codon C-G-G-A/G.
- The apparent "CGG→Gly" from Codetta/KACI-type methods is the in-frame shadow of this.
- No bacterial genome-wide quadruplet/frameshift code has been reported. Closest precedents: Perkinsus mitochondria
  (unused-codon-triggered frameshifts) and Euplotes (stop-codon +1 PRF).

**Critical tests pending (VM down):**
1. All CGG sites — next-nt distribution (expect A/G), and frame 0 versus +1 homology and stop distance, CAG-508
   versus control.
2. Cross-genus conservation of individual shift sites; sequencing errors would not be shared.
3. Are other single-nt insertions elevated in CAG-508 assemblies (general error) or only C-GGR?
4. Correlation with presence of the tRNA-Gly G37 variant.

## Caveats so far
- Mostly MAGs; no isolates of CAG-508 are known to us.
- The per-genome signal is weak because CGG is very rare; inference relies on pooling.
- The decoding mechanism is unknown. Ambiguous decoding or partial "unassigned" status cannot be excluded.
- The tRNA survey uses ARAGORN anticodon labels plus tRNAscan for a subset; isotype identity of the Gly tRNAs is
  not resolved.
