#!/usr/bin/env python
"""Build the 4-6 page executive summary PDF (reportlab Platypus)."""
import os, matplotlib
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, KeepTogether, PageBreak)
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
fd = os.path.join(os.path.dirname(matplotlib.__file__), 'mpl-data', 'fonts', 'ttf')
pdfmetrics.registerFont(TTFont('DV', os.path.join(fd, 'DejaVuSans.ttf')))
pdfmetrics.registerFont(TTFont('DVB', os.path.join(fd, 'DejaVuSans-Bold.ttf')))
pdfmetrics.registerFont(TTFont('DVI', os.path.join(fd, 'DejaVuSans-Oblique.ttf')))
pdfmetrics.registerFont(TTFont('DVM', os.path.join(fd, 'DejaVuSansMono.ttf')))
from reportlab.pdfbase.pdfmetrics import registerFontFamily
registerFontFamily('DV', normal='DV', bold='DVB', italic='DVI', boldItalic='DVB')
INK = colors.HexColor('#0b0b0b'); INK2 = colors.HexColor('#52514e'); BLUE = colors.HexColor('#2a78d6'); RULE = colors.HexColor('#d9d7d0')
S = dict(
    title=ParagraphStyle('t', fontName='DVB', fontSize=14.5, leading=18, textColor=INK, spaceAfter=3),
    sub=ParagraphStyle('s', fontName='DV', fontSize=8.5, leading=11, textColor=INK2, spaceAfter=6),
    h=ParagraphStyle('h', fontName='DVB', fontSize=10.5, leading=13, textColor=BLUE, spaceBefore=7, spaceAfter=3),
    b=ParagraphStyle('b', fontName='DV', fontSize=8.4, leading=11.2, textColor=INK, spaceAfter=3),
    bl=ParagraphStyle('bl', fontName='DV', fontSize=8.4, leading=11.2, textColor=INK, leftIndent=10, bulletIndent=2, spaceAfter=1.5),
    cap=ParagraphStyle('c', fontName='DV', fontSize=7.3, leading=9.4, textColor=INK2, spaceAfter=6),
    tc=ParagraphStyle('tc', fontName='DV', fontSize=7.2, leading=9, textColor=INK),
    th=ParagraphStyle('th', fontName='DVB', fontSize=7.2, leading=9, textColor=INK),
)
P = lambda t, s='b': Paragraph(t, S[s])
def bullets(items): return [Paragraph(i, S['bl'], bulletText='•') for i in items]
def table(rows, widths, head=True):
    data = [[Paragraph(str(c), S['th'] if (head and r == 0) else S['tc']) for c in row] for r, row in enumerate(rows)]
    t = Table(data, colWidths=widths, repeatRows=1 if head else 0)
    t.setStyle(TableStyle([('LINEBELOW', (0, 0), (-1, 0), 0.8, INK2), ('LINEBELOW', (0, 1), (-1, -1), 0.3, RULE),
                           ('VALIGN', (0, 0), (-1, -1), 'TOP'), ('TOPPADDING', (0, 0), (-1, -1), 2.2), ('BOTTOMPADDING', (0, 0), (-1, -1), 2.2),
                           ('LEFTPADDING', (0, 0), (-1, -1), 3), ('RIGHTPADDING', (0, 0), (-1, -1), 3)]))
    return t
def fig(path, w_cm, cap):
    from PIL import Image as PI
    im = PI.open(path); r = im.size[1] / im.size[0]
    return KeepTogether([Image(path, width=w_cm * cm, height=w_cm * cm * r), P(cap, 'cap')])
F = 'figs/'
st = []
st += [P('A conserved CGG-triggered +1 frameshift code in the gut bacterial order TANB77, coinciding with loss of every CGG-reading tRNA and gain of a tRNA<super>Gly</super> variant', 'title'),
       P('Executive summary · computational discovery from public genomes (GTDB r226; 143,614 species) · 20 Nov 2026 · '
         'Full report and code: <font name="DVM">github.com/rithinnagaraj/tanb77-cgg-frameshift</font>', 'sub')]
st += [P('1. Summary', 'h'),
       P('We screened all 143,614 GTDB r226 bacterial and archaeal species for non-standard genetic codes. The screen recovered every known bacterial code '
         'and found one new anomaly: in the human/animal gut family <b>CAG-508</b> (order <b>TANB77</b>, Clostridia; 830 species; '
         'a taxon associated with cancer-immunotherapy response), the arginine codon <b>CGG</b> appears to code for glycine. '
         'Close inspection showed this is not a reassignment. CGG is a conserved <b>+1 frameshift site</b>: the protein stays homologous only if the C is skipped '
         'and the following GGA/GGU is read as Gly. At 74.5% of 1,284 CGG sites the +1 frame matches homologs better (median +57.9 bits); '
         'control codons in the same genomes and CGG in other Clostridia do so at about 1%. '
         'The sites are inherited, not assembly artefacts. They are shared across genera (one site in Mfd occurs in 84 genomes from 13 genera) and present in '
         '<b>23 complete genomes from 6 independent projects, 19 of them PacBio HiFi</b>. '
         'The phenomenon coincides with the loss of every CGG-reading tRNA (anticodons CCG/UCG) and the gain of a single-copy '
         '<b>tRNA<super>Gly</super>(UCC) variant with G37</b>, whose loop could pair a C-GGA quadruplet. The variant is present in all 23 complete genomes, '
         '65–96% of TANB77 MAGs, and 0 of 305 control genomes. '
         'Automated code-inference tools mistake these sites for a "CGG→Gly" reassignment (Melnykov\'s public KACI table holds 863 such undescribed calls). '
         'Standard gene annotations of TANB77 are systematically truncated at CGG. '
         'To our knowledge this is the first clade-wide, tRNA-linked use of a sense codon as a conserved frameshift signal in a prokaryote.'),
       fig(F + 'fig0_model.png', 15.5, '<b>Figure 1. Model.</b> Prodigal/GTDB read CGG in frame (orange), which loses homology and hits an early stop. '
           'The proposed reading (blue) consumes C-GGA as one Gly and continues in +1. The sequence is TIGR02386 (RNA polymerase β′) from CAG-508.')]
st += [P('2. Specific claims', 'h'),
       table([['#', 'Claim', 'Key evidence', 'Confidence'],
              ['C1', 'CGG in conserved CAG-508 genes is predominantly a +1 frameshift site', '74.5% of 1,284 sites prefer +1 (controls 0.7–1.3%); position-matched gene tails: 82% vs 8.6%; interior native sites 91.5% (n=117)', 'High'],
              ['C2', 'The shift sites are real, inherited sequence, not sequencing or assembly error', 'Present in 23 complete genomes (19 HiFi, 6 projects); 73/281 sites shared across ≥2 genera; C-before-GGN is 16× more frequent than A and 110× more than U', 'High'],
              ['C3', 'It coincides with loss of all CGG-reading tRNAs in TANB77 (one evolutionary event)', 'Reader in 0.5% of CAG-508 and 0/23 complete genomes; outgroup RGIG9151 keeps a reader (5/5) and does not shift', 'Moderate–high (single contrast)'],
              ['C4', 'A conserved tRNA<super>Gly</super>(UCC) G37 variant co-occurs and plausibly mediates C-GGA decoding', 'Variant in 23/23 complete genomes (single copy), 0/300 other Clostridia; CGGA over-represented vs GGA usage (72% vs 59%, z≈9.4)', 'Correlative; mechanism hypothesised'],
              ['C5', 'Standard annotations truncate many TANB77 genes; the "CGG→Gly" calls of code-inference tools are a frameshift shadow', '42.5% of CAG-508 CGGs sit in the last 15 codons of genes (vs 5.3%); 66% have a +1 ORF continuing (vs 13–14%); KACI and our screen both call G', 'High'],
              ['C6', 'A minority of CGG is read in frame', 'For example, 43 genomes carry CGG at an invariant Arg column of purB without a shift (~4% of native sites)', 'Moderate']],
             [0.8 * cm, 5.3 * cm, 7.5 * cm, 2.6 * cm])]
st += [fig(F + 'fig1_frame_test.png', 15.5, '<b>Figure 2.</b> (A) Genome-level frame test. For every CGG and a sample of control GGN codons in marker genes, the CDS was located in the genome, and '
           '40 aa upstream plus 70 aa downstream were scored against the marker HMM in frame 0 versus frame +1. (B) The base after CGG. After the +1 shift, CGG-A/U becomes a GGA/GGU glycine codon '
           '(CAG-508\'s own GGN third-position usage is A 59%, U 32%).'),
       fig(F + 'fig2_alignment.png', 15.5, '<b>Figure 3.</b> Two of the most widely shared shift sites. The extra C (bracketed) sits immediately before GGA; reading +1 restores the conserved sequence.')]
st += [P('3. Novelty: how it was checked', 'h')] + bullets([
    '<b>Literature</b> (Exa web search through Nov 2026). Genetic-code screens: Shulgina & Eddy 2021; Fakih 2026; Melnykov 2026; Kim 2026 bioRxiv; gTranslate and Parks 2026; Kivenson 2025 (Pyl code). '
    'Frameshift and quadruplet decoding: DeBenedictis 2022 notes the "absence of naturally occurring quadruplet codon codes"; Euplotes; Perkinsus mitochondria; SufB2/SufD; FSscan 2009; Antonov 2013. '
    'TANB77: Kim et al. 2024. None reports CGG-triggered or clade-wide frameshifting in TANB77 or in any bacterium.',
    '<b>Data-level checks against prior screens.</b> Codetta 2021 supplementary files: 68 CAG-508 genomes were called CGG \'?\' or \'R\', and there is no candidate. '
    'Melnykov\'s full KACI table (3.1M assemblies, Zenodo 20837797): 863 CAG-508 assemblies have CGG=G, which is 98.6% of all GTDB-mapped CGG→Gly calls. These were automated and never described; his paper defers '
    'Clostridial CGG decoding to "a separate investigation". We show the signal is a frameshift shadow.',
    '<b>Precedents that narrow the claim.</b> Isolated +1 shifts at rare or uncharged codons are known: E. coli CGG-C/AGG-AGG sites (Spanjaard & van Duin 1988; FSscan), recoded E. coli (Ma 2018), '
    'and yeast near-cognate codons (Sundararajan 1999). Eukaryotic organelle and ciliate frameshift codes are known (Perkinsus, Euplotes). Our novelty is a clade-wide, conserved, tRNA-linked frameshift code in bacteria.',
    'Two independent adversarial-review agents searched for prior art and found none.'])
st += [fig(F + 'fig5_longread_trna.png', 15.5, '<b>Figure 4.</b> (A) Replication in 23 complete CAG-508 genomes (PacBio HiFi, Nanopore hybrid, Illumina; 6 BioProjects). '
           'CGG gene-tail sites continue in +1 far more often than control AGG/GGG sites. (B) TANB77 lineages that lost the CGG-reading tRNA carry a tRNA<super>Gly</super>(UCC) G37 variant; '
           'outgroups have the reader and no variant (MAG tRNA recovery is incomplete; complete genomes: 23/23 variant, 0/23 reader).')]
st += [P('4. Testable predictions', 'h')] + bullets([
    '<b>Proteomics / ribosome profiling.</b> In enriched or cultured TANB77, peptides will span C-GGA junctions with Gly at the shift when searched against shift-corrected proteins, and the annotated in-frame tails will be absent. '
    'Ribosome profiling will show pausing at CGG. Public metaproteomes cannot test this yet: MetaPep holds ≤5 peptides per CAG-508 species.',
    '<b>Reporter / tRNA test.</b> A CAG-508 C-GGA site in E. coli (which has tRNA<super>Arg</super>CCG) will be read in frame as Arg and terminate. '
    'Expressing the TANB77 tRNA<super>Gly</super>(UCC) G37 variant, especially in a ΔtRNA-CCG background, will induce +1 decoding with Gly incorporation; the A37 paralog will not.',
    '<b>Genotype rule.</b> Any TANB77 genome that re-acquires a CGG-reading tRNA should lose the shift signature; RGIG9151 fits this already. Other lineages that lose all CGG readers and gain an equivalent tRNA should evolve the same signature. '
    'A Peptostreptococcales MAG (GCA_937907205, 55% shift fraction) is a candidate.',
    '<b>Context.</b> Shift sites will be biased to C-GGA over C-GGU relative to GGN usage, as the UCC anticodon of the variant implies; observed 72% vs 59%.'])
st += [P('5. Controls for the method: known results and false flags', 'h'),
       table([['Control', 'Expectation', 'Result'],
              ['Known bacterial codes (4, 25, 37; RFN20 CGG→Trp; Enterosoma AGG→Met)', 'Recovered', 'Recovered (427 UGA→W, 204 UGA→G, 56 CGA→W, 4 CGG→W, 1 AGG→M genome calls)'],
              ['Control GGN codons in the same CAG-508 genomes', 'No +1 preference', '1.3% (78,544 sites); gene-tail positions 8.6%'],
              ['CGG in other Clostridia', 'No +1 preference', '0.7% (54,478 sites)'],
              ['Insertion artefact: X-GG at Gly columns per GGN (X = A, C, U)', 'Similar for all X if it is an error', 'CAG-508: C 5.0e-4, A 3.2e-5, U 4.5e-6; A and U are at the other-Clostridia level'],
              ['All 13 GTDB families with rare CGG (<5 per 10k codons)', 'CGG at Arg sites unless reassigned', 'Arg in 97–100% except the known Trp clades and CAG-508 (Fig. 5B)'],
              ['Calibrated null over 812 rare codon×family combinations', 'Rarely ≥0.85 non-standard', 'p=0.0012; the only exceedance is an Ile/Leu similarity artefact'],
              ['False flags caught', '', 'Pooled-only ATG/AGY "reassignments" with no supporting genome; GC-rich contaminant gene copies carrying Arg-CGG (98% of Arg-site CGGs); Codetta 2021 \'R\' calls driven by contamination; our own initial "CGG→Gly" call']],
             [5.2 * cm, 3.6 * cm, 7.4 * cm])]
st += [fig(F + 'fig3_trna_link.png', 15.5, '<b>Figure 5A.</b> Frameshift signature by TANB77 family (genome-wide, all CDSs). Only families that lost every CGG-reading tRNA show it.'),
       fig(F + 'fig4_specificity.png', 15.5, '<b>Figure 5B.</b> (A) CGG at Gly columns per source codon, against the CGN controls. (B) In every family where CGG is rare, CGG sits at Arg columns; '
           'the exceptions are the known CGG→Trp clades and CAG-508, where the "Gly" columns are the shift sites.')]
st += [P('6. Verification ladder', 'h'),
       table([['Rung', 'Test', 'Outcome'],
              ['1', 'Marker screen; pooled per clade; one genome per genus (44 genera, 5 resamples); CAG-269 left out', 'Anomaly robust: P=1 in every resample'],
              ['2', 'Independent whole-proteome Codetta re-implementation (all Pfam, markers excluded)', 'Same CGG signal in all 5 genus representatives'],
              ['3', 'Independent method and author (KACI, 3.1M assemblies)', '863 CAG-508 assemblies; 98.6% of all GTDB-mapped calls'],
              ['4', 'Contamination split by gene GC', 'Arg-site CGGs sit on foreign GC-rich copies (98%)'],
              ['5', 'Frame test on genomes, with position-matched controls', '+1 restores homology at CGG only (74.5% vs about 1%)'],
              ['6', 'Corrected insertion comparator (AGG/UGG)', 'Specific to C (16–110×)'],
              ['7', 'Cross-genus site sharing', '73 sites shared; one site in 13 genera'],
              ['8', 'Long-read / independent-platform replication (23 complete genomes)', 'Signature in all 23; 19 PacBio HiFi'],
              ['9', 'tRNA genotype–phenotype; G37 variant survey (1,228 genomes)', 'Reader loss and variant gain coincide with the shift; outgroups negative'],
              ['10', 'Two adversarial reviews, with re-analysis after each', 'Comparator, null and wording corrected; long-read test added'],
              ['—', 'Not done: translation-level evidence (proteomics, ribosome profiling)', 'Public data are insufficient; prediction stated']],
             [1.1 * cm, 8.4 * cm, 6.7 * cm])]
st += [P('7. Limitations', 'h')] + bullets([
    'No experimental or translation-level evidence. The amino acid at the shift (Gly) is inferred from the +1 codon, not observed. The mechanism (quadruplet pairing by the G37 tRNA versus P-site slippage) is a hypothesis.',
    'TANB77 has no isolates. The complete genomes are metagenome-derived, though largely HiFi and from 6 projects.',
    'tRNA absence is called from assemblies (ARAGORN). It is reliable in complete genomes and incomplete in MAGs, and the tRNA loss is a single evolutionary contrast.',
    'About 25% of CGG sites do not prefer +1. These are partly contaminant copies, and partly in-frame sites such as purB. Shift efficiency, regulation and fitness cost are unknown.',
    'The frame test uses marker genes plus a genome-wide ORF heuristic. A frameshift-aware census of all proteins is still to be done.'])
st += [P('8. Implications', 'h')] + bullets([
    '<b>Genetic code evolution.</b> When codon capture empties a codon of its tRNA, the codon can be kept as a conserved frameshift signal, apparently with a co-evolved tRNA, rather than reassigned to a new amino acid. '
    'This is a bacterial counterpart to Euplotes and Perkinsus.',
    '<b>Methods.</b> Codetta-style and k-mer code inference can report spurious sense reassignments at frameshift codons. Frame-aware checks should accompany any reassignment claim in clades that lack the cognate tRNA, '
    'for example the other Clostridial CGG→Gln/Trp/Ala calls.',
    '<b>Microbiome science.</b> Gene catalogues (GTDB, UHGG, MGnify) truncate dozens of genes per TANB77 genome. TANB77 is linked to immunotherapy response, so its functional annotation and metaproteomics need shift-aware gene models.',
    '<b>Synthetic biology.</b> This is a naturally evolved C-GGA decoding context and tRNA, and a candidate component for quadruplet-codon engineering.'])
st += [P('9. Method in brief and code', 'h'),
       P('<b>Screen.</b> GTDB r226 bac120/ar53 marker CDS for all species representatives were aligned to Pfam/TIGRFAM HMMs with pyhmmer. Codons were mapped to match columns, and Codetta-normalised '
         'likelihoods were computed per genome and pooled by clade. <b>Confirmation.</b> Whole-proteome Pfam hmmsearch with column→codon mapping; the frame test on genomes; a genome-wide CDS-tail and +1-ORF census; '
         'ARAGORN/tRNAscan-SE tRNA profiles; 23 complete NCBI genomes gene-called with pyrodigal. <b>Novelty data</b>: Codetta 2021 supplementary files and the Melnykov KACI table. '
         '<b>Code &amp; data</b>: <font name="DVM">https://github.com/rithinnagaraj/tanb77-cgg-frameshift</font> (code/: gc_screen, gc_infer, gc_genome, gc_checks, gc_sites, gc_frames_genome, gc_genomewide, gc_profile, gc_glyvar, gc_comparator, make_figs, build_pdf). '
         'All inputs are public.')]
st += [P('10. Key references', 'h'),
       P('Shulgina & Eddy, eLife 10:e71402 (2021) · Melnykov, PNAS 123:e2610659123 (2026); Zenodo 10.5281/zenodo.20837797 · Fakih et al., NAR 54:gkag327 (2026) · '
         'Kim et al., bioRxiv 2026.11.10.764470 (2026) · Parks et al., Microb Genom 12:001767 (2026) · Chaumeil et al., PeerJ (2026, gTranslate) · Kivenson et al., Science 390:eadu2404 (2025) · '
         'DeBenedictis et al., eLife (2022) · Gornik et al., GBE (2022, Perkinsus) · Wang et al., BMC Genomics 26:676 (2025, Euplotes) · Sundararajan et al., Mol Cell (1999) · Spanjaard & van Duin, PNAS (1988) · '
         'Ma et al., eLife (2018) · Antonov et al., NAR (2013) · Oba et al., PNAS 88:921 (1991) · Kim et al., Nat Commun (2024, TANB77) · Sun et al., MCP 23:100763 (2024, MetaPep) · Parks et al., NAR 54:D743 (2026, GTDB r226).', 'cap')]
st += [PageBreak(), P('Appendix A (optional). Adversarial review summary', 'h')] + bullets([
    '<b>Plan review (riboswitch project).</b> Required a COMPASS-atlas novelty gate. The gate found only known <i>pan</i>-motif variants and BirA operators, so we pivoted to a genetic-code screen.',
    '<b>Review 1 ("CGG→Gly").</b> Estimated about 40% for Gly decoding and flagged the denominator, pseudo-replication, the Gly/Ala question and the Melnykov data. The re-analysis passed every check, and the alignment inspection prompted '
    'the frameshift reinterpretation.',
    '<b>Review 2 (frameshift).</b> Estimated P(C1+C2) about 0.75 and P(C3 causal) about 0.6. It flagged the wrong insertion comparator, the next-base null, unmatched controls, single-platform data and overclaiming wording. '
    'All were addressed: AGG/UGG comparator, the CAG-508 GGN null, position-matched controls, 23 long-read genomes, the tRNA variant and recalibrated wording. '
    'Still open: read-level SRA mapping and a frameshift-aware census of all proteins.'])
st += [P('Appendix B (optional). Human vs AI contribution', 'h'),
       P('<b>Human (user):</b> set the goal and the novelty bar, provisioned and scaled the GCP VM, rebooted it after memory hangs, approved the download of public data, and relayed the requirement that a negative result is insufficient, '
         'plus two requested checks (interior/terminal split, proteomics feasibility). <b>AI (Claude):</b> chose questions, ran literature and novelty checks, made the pivot decisions, designed and wrote all code, ran the analyses, '
         'interpreted results (including retracting its own "CGG→Gly" interpretation), ran independent reviewer agents, and wrote this report.')]
doc = SimpleDocTemplate('report/executive_summary.pdf', pagesize=A4, leftMargin=1.6 * cm, rightMargin=1.6 * cm, topMargin=1.4 * cm, bottomMargin=1.4 * cm,
                        title='CGG +1 frameshift code in TANB77 — executive summary', author='Claude (with the user)')
def foot(c, d):
    c.saveState(); c.setFont('DV', 7); c.setFillColor(INK2); c.drawRightString(A4[0] - 1.6 * cm, 0.8 * cm, f'{d.page}'); c.restoreState()
doc.build(st, onFirstPage=foot, onLaterPages=foot)
from pypdf import PdfReader
print('pages', len(PdfReader('report/executive_summary.pdf').pages))
