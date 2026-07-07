# HAR annotator

A reproducible pipeline that filters **Human Accelerated Regions (HARs)** to
constrained, neurodevelopment-adjacent, disease-overlapping elements, assembles
a per-element evidence table, and ranks candidates with a **transparent
additive score** — then interprets the top neurodevelopmental candidates with
an antagonistic-pleiotropy framing.

Everything is hg38. Every data pull is cached, hashed, and logged to
`data/manifest.csv`, so the analysis replays end-to-end from source.

```
python run_pipeline.py            # reproduces the shipped shortlist exactly
```

---

## What it does — the funnel

| Step | Filter | HARs |
|------|--------|------|
| 0 | All HARs (Cui et al. 2025, *Nature*, hg38) | 3,257 |
| A | Mammalian constraint — mean 241-way phyloP > 1.0 | 2,757 |
| B | Within 1 Mb of a neurodevelopmental-disorder gene TSS | 1,718 |
| C | Within 25 kb of a genome-wide-significant neuropsychiatric GWAS lead SNP | **363** |

The 363 candidates are then scored on six axes and ranked. Top candidate:
**HAR_575 → *POC1B*** (PLAC-seq-linked, educational-attainment GWAS overlap,
embryonic-brain DNase-active). Classic neurodevelopmental / human-evolution
genes recovered in the shortlist include *ZEB2, TCF4, MEF2C, PHOX2B, TCF20,
ZSWIM6, FOXP2, SOX5, NR4A2, MITF*.

---

## The transparent score

The total is a **weighted sum of six normalized (0–1) component scores** — no
black box. Every `score_<c>` and its weighted `contrib_<c>` is retained as a
column, so `total_score` decomposes exactly.

```
total_score = Σ  WEIGHTS[c] · score_c        (c in the six components below)
```

| Component | Default weight | Evidence |
|-----------|:-:|----------|
| `constraint`   | 0.20 | mean 241-way phyloP, min-max scaled |
| `gene`         | 0.25 | DDG2P confidence tier × TSS-distance decay, +0.15 if PLAC-seq-linked |
| `disease`      | 0.25 | −log10(GWAS p) × proximity to lead SNP |
| `brain`        | 0.15 | ENCODE embryonic-brain DNase peak overlap (+ PLAC-seq ATAC support) |
| `acceleration` | 0.10 | HAR width (log-scaled) — **proxy**, see caveats |
| `motif`        | 0.05 | reserved (0 unless a motif-disruption annotation is supplied) |

Weights are module constants (`har_annotator/score.py`) and overridable per run:
`python run_pipeline.py --weights gene=0.35,disease=0.30`.

---

## Install & run

```bash
# conda env (Python 3.11 — pyranges pins <3.13)
conda create -n hars -c bioconda -c conda-forge python=3.11 \
    pandas numpy scipy matplotlib seaborn requests \
    pyranges=0.1.4 pybigwig biopython pyarrow openpyxl
conda activate hars

python run_pipeline.py                     # defaults, reproduces shipped outputs
python run_pipeline.py --gwas-window 50000 # sensitivity: wider disease window
python run_pipeline.py --min-phylop 1.5 --gene-window 500000
python run_pipeline.py --weights disease=0.35,gene=0.20 --outdir alt_weights/
```

Key parameters (all have defaults that reproduce the shipped shortlist):

| Flag | Default | Meaning |
|------|:-:|---------|
| `--min-phylop`   | 1.0 | Step A constraint threshold (mean phyloP) |
| `--gene-window`  | 1_000_000 | Step B max HAR–TSS distance (bp) |
| `--gwas-window`  | 25_000 | Step C max HAR–lead-SNP distance (bp) |
| `--gwas-pval`    | 5e-8 | GWAS genome-wide significance cutoff |
| `--weights`      | — | override score weights, `k=v,k=v` |
| `--sfari`        | — | optional SFARI gene CSV, unioned into the neurodev set |
| `--force-download` | off | re-fetch all sources instead of using the cache |

GWAS-window sensitivity (candidate count): 10 kb → 196, **25 kb → 363**,
50 kb → 526, 100 kb → 759. 25 kb (LD-block scale) is the shipped default.

---

## Data sources

All fetched, versioned, and hashed in `data/manifest.csv` (URL + access date +
SHA-256 + size). See `phase0_sources.md` for the full provenance / de-risking
notes.

| Source | Used for | manifest row | builder / consumer |
|--------|----------|:-:|:-:|
| **HARs & PLAC-seq** — Cui et al. 2025, *Nature* | 3,257 HARs, hg38-native; neuronal PLAC-seq HAR→gene links | `hars_bed`, `hars_meta`, `cui2025_supp4` | `data_io.py`, `references.build_plac_links` |
| **Constraint** — Zoonomia 241-way phyloP | per-HAR mean/max phyloP | *(remote query, not cached)* | `filters.annotate_phylop` |
| **Neurodev genes** — DDG2P (Gene2Phenotype) | 2,524 genes → hg38 coords + confidence tier | `ddg2p`, `refgene_hg38` | `references.build_neurodev` |
| **GWAS** — EBI GWAS Catalog (ontology-annotated) | 22,489 genome-wide-sig neuropsychiatric lead SNPs | `gwas_assoc` | `references.build_gwas_loci` |
| **Developing brain** — ENCODE embryonic-brain DNase-seq | 165,568 peaks; regulatory-activity axis | `fetal_brain_dnase` | `evidence.annotate_brain_dnase` |
| **Gene coords** — UCSC refGene hg38 | gene symbol → strand-aware TSS | `refgene_hg38` | `references.build_genes` |
| **TF motifs** | *(reserved — not wired in; `motif` score = 0)* | — | — |

Every raw file's exact download URL, SHA-256, access date, and size is in
`data/manifest.csv`; `phase0_sources.md` has the full de-risking notes. The
raw files themselves are **not committed** (see `.gitignore`) — `download.py`
re-fetches and hash-verifies them from the URLs below.

### Primary references & direct links

- **HARs & PLAC-seq map — Cui et al. (2025), *Nature* 640:991–999.**
  "Comparative characterization of human accelerated regions in neurons,"
  doi:[10.1038/s41586-025-08622-x](https://doi.org/10.1038/s41586-025-08622-x).
  Supplies both the 3,257 HAR coordinates **and** the neuronal PLAC-seq HAR→gene
  map — they are the same table, not separate downloads. That table is the
  paper's **Supplementary Table 2** ("HAR and their chimpanzee ortholog
  interacting genes in neurons"); the map used here is its "Table-b: Genes
  interacting with HARs in human" sheet (1,303 HARs, 1,719 genes,
  genome-wide). **Note on numbering:** the download file is named
  `41586_2025_8622_MOESM4_ESM.xlsx` (Nature's *MOESM4* file id), but the table
  it contains is Supplementary Table **2** — the MOESM file number and the
  supplementary-table number do not correspond. This is *not* Supplementary
  Table 3, which is the separate CRISPRi-prioritized shortlist of ~20 HARs
  selected for functional validation. Fetched via the GitHub mirror
  [`athenamarou/HAR-TFBS-Project`](https://github.com/athenamarou/HAR-TFBS-Project)
  (`data/hars_hg38.bed`, `data/hars_hg38.tsv`,
  `data/supplementary/41586_2025_8622_MOESM4_ESM.xlsx`).

  **Method note:** the HAR→gene links come from **PLAC-seq** (proximity
  ligation-assisted ChIP-seq, anchored on the H3K4me3 active-promoter mark) in
  human and chimpanzee iPSC-derived neurons — a ChIP-anchored chromatin-contact
  assay, *not* generic Hi-C. Earlier drafts of this project called it "Hi-C";
  that was imprecise and has been corrected throughout, including the column
  names (`plac_gene`, `plac_type`, `plac_atac`, assignment `plac_linked`).
- **Neurodevelopmental gene list — DDG2P, from Gene2Phenotype (G2P).**
  Download: `https://ftp.ebi.ac.uk/pub/databases/gene2phenotype/G2P_data_downloads/2026_06_28/DDG2P_2026-06-28.csv.gz` (EBI).
  Portal: [www.ebi.ac.uk/gene2phenotype](https://www.ebi.ac.uk/gene2phenotype).
  Primary paper: Thormann *et al.* (2024), *Genome Medicine* — "Curating genomic
  disease-gene relationships with Gene2Phenotype (G2P)",
  doi:[10.1186/s13073-024-01398-1](https://doi.org/10.1186/s13073-024-01398-1)
  ([PMC11539801](https://pmc.ncbi.nlm.nih.gov/articles/PMC11539801/)); earlier
  G2P/VEP tool paper: Thormann *et al.* (2019), *Nature Communications*,
  doi:[10.1038/s41467-019-10016-3](https://doi.org/10.1038/s41467-019-10016-3).
- **Developing-brain open chromatin — ENCODE DNase-seq `ENCFF660HML`.**
  Peak file: [encodeproject.org/files/ENCFF660HML](https://www.encodeproject.org/files/ENCFF660HML/)
  (download: `@@download/ENCFF660HML.bed.gz`); parent experiment
  [ENCSR420RWU](https://www.encodeproject.org/experiments/ENCSR420RWU/)
  (brain, male embryo 105 days). Portal paper: Davis *et al.* (2018),
  *Nucleic Acids Research* — "The Encyclopedia of DNA Elements (ENCODE): data
  portal update", doi:[10.1093/nar/gkx1081](https://doi.org/10.1093/nar/gkx1081);
  consortium: The ENCODE Project Consortium (2012), *Nature*,
  doi:[10.1038/nature11247](https://doi.org/10.1038/nature11247).
- **Neuropsychiatric GWAS — NHGRI-EBI GWAS Catalog** (ontology-annotated
  associations, latest release).
  Download: `https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip`;
  portal [www.ebi.ac.uk/gwas](https://www.ebi.ac.uk/gwas).
  Papers: Cerezo *et al.* (2025), *Nucleic Acids Research* 53(D1):D998–D1005,
  doi:[10.1093/nar/gkae1070](https://doi.org/10.1093/nar/gkae1070) (current
  release); Sollis *et al.* (2023), *NAR* 51(D1):D977–D985,
  doi:[10.1093/nar/gkac1010](https://doi.org/10.1093/nar/gkac1010).
  **"Neuropsychiatric" is a keyword filter, not a Catalog category.**
  `references.build_gwas_loci` selects rows whose `MAPPED_TRAIT`/`DISEASE/TRAIT`
  match `references._NEURO_PAT` (schizophren·autis·bipolar·depress·neurotic·
  cognit·educational·intellig·attention-deficit/adhd·alzheimer·neurodevelop·
  intellectual-disab·epileps·tourette·obsessive) and keeps genome-wide-
  significant hits (p < 5×10⁻⁸) → 22,489 lead SNPs. The term list is a module
  constant and is meant to be edited. **Trait-mix caveat:** the largest
  contributors are *educational attainment*, *intelligence*, and *cognitive-
  function* GWAS — behavioural/cognitive proxies with known social-environmental
  confounding. "Overlaps an educational-attainment locus" is therefore a weaker
  biological claim than "overlaps a schizophrenia locus"; the shipped
  `mapped_trait` column lets you see which applies to each candidate.
- **Constraint — Zoonomia 241-way phyloP** (Cactus 241-mammal alignment),
  queried remotely from the UCSC bigWig (`hg38 cactus241way`); 9.6 GB, never
  downloaded. The phyloP scores themselves: Christmas *et al.* (2023),
  *Science* 380:eabn3943, "Evolutionary constraint and innovation across
  hundreds of placental mammals,"
  doi:[10.1126/science.abn3943](https://doi.org/10.1126/science.abn3943).
  Disease-mapping utility of base-pair mammalian constraint: Sullivan *et al.*
  (2023), *Science* 380:eabn2937,
  doi:[10.1126/science.abn2937](https://doi.org/10.1126/science.abn2937).
- **Gene coordinates — UCSC refGene (hg38).**
  `https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz`.

**SFARI substitution note.** SFARI Gene was unreachable from the build
environment (`gene-archive.sfari.org` blocked), so the neurodevelopmental gene
set is **DDG2P only** — the shipped 363 candidates reflect that. The pipeline
unions any SFARI symbols supplied via `--sfari`, so a SFARI export can be
dropped in without code changes (this will change the candidate count and
ranking).

**TF-motif note.** The `motif` score component (weight 0.05) is **reserved and
contributes 0** — no motif database is downloaded or wired in. Scoring whether
each HAR's human-specific substitutions create/disrupt a TF binding site
(e.g. against JASPAR / HOCOMOCO PWMs) is flagged as future work, alongside the
true human/chimp substitution count.

---

## Package layout

```
har_annotator/
  download.py    cached + hashed + manifest-logged fetches
  data_io.py     load HAR base table + Cui PLAC-seq sheet from the supplement
  references.py  build derived tables: genes, neurodev (DDG2P), GWAS loci, PLAC-seq links
  filters.py     Phase-1 funnel: annotate_phylop, filter_constrained,
                 assign_nearest_gene, annotate_gwas, keep_gwas
  evidence.py    Phase-2 per-element evidence spine (assemble)
  score.py       Phase-3 transparent additive score (compute_scores, WEIGHTS)
run_pipeline.py  end-to-end driver (parameterized)
```

Outputs: `funnel_counts.csv`, `candidate_hars.parquet`, `har_evidence.parquet`
(+ `har_evidence_schema.json`), `har_shortlist_ranked.csv`/`.parquet`,
`har_shortlist_top50.csv`. Interpretation: `top_candidates_interpretation.md`.
Figures: `fig_funnel.png`, `fig_score_overview.png`,
`fig_top_candidate_detail.png`.

---

## What "mammalian constraint" means here

**phyloP** (phylogenetic P-value) is a per-base evolutionary-rate score. For
each position in the human genome, it compares the substitution rate observed
across an alignment of 241 placental mammals against the rate expected under
neutral drift. Positive = **slower than neutral** = the base has been conserved
= *constrained* (purifying selection removed mutations there). Negative = faster
than neutral (accelerated); ~0 = evolving neutrally. The scores come from the
Zoonomia Cactus 241-mammal alignment (Christmas et al. 2023). We take the
**mean phyloP over each HAR's span** as its constraint value and keep HARs with
mean > 1.0 (Step A: 3,257 → 2,757).

Why constraint matters for a *human-accelerated* region — the apparent paradox
is the point. A HAR is a locus that was **deeply conserved across mammals**
(high ancestral constraint) yet **changed rapidly on the human branch**. The
mammalian-constraint filter enriches for elements that were doing something
selectively important for a long time; the human-specific change is then more
likely to be functionally consequential than a change in a locus that was
never constrained. Constraint (mammalian, phyloP) and acceleration
(human-branch) are **different axes** — this pipeline scores constraint
directly and, as noted in the caveats, only *proxies* acceleration.

**"Sharpened disease mapping" — what's actually established vs. what we
verified here.** The strong, external result is that base-pair mammalian
constraint from Zoonomia enriches for trait/disease heritability and improves
prioritization of causal variants — shown genome-wide in Sullivan et al. (2023,
*Science*, doi:10.1126/science.abn2937), not something this pipeline proves. To
avoid overclaiming, the enrichment was checked **inside this dataset**: of all
3,257 HARs, constrained ones (phyloP > 1) overlap a neuropsychiatric GWAS locus
(±25 kb) at **21.0%** vs **15.4%** for unconstrained — a **1.37× enrichment**
(Fisher odds 1.46, p = 3.6×10⁻³). So constraint does concentrate disease
overlap in our HAR set, but **modestly**; it is a sensible prioritization axis,
not a strong causal sieve on its own. The pipeline's power comes from
*combining* it with gene, disease, and brain-activity evidence.

---

## Interpretation & caveats

Read `top_candidates_interpretation.md` for the top-candidate biology. Three
caveats apply throughout and are load-bearing:

1. **Overlap ≠ causation.** "Overlaps a GWAS locus" means within 25 kb of a
   lead SNP — not that the HAR is the causal element or that its
   human-specific changes drive the association.
2. **Nearest gene ≠ target gene.** `nearest_tss` assignments are a proximity
   heuristic; `plac_linked` assignments rest on a neuronal PLAC-seq interaction (stronger,
   still correlative).
3. **Acceleration is a proxy here.** The Cui table ships no per-element
   acceleration statistic we could thread through, so HAR width (log-scaled)
   stands in and mammalian phyloP carries the constraint axis. A true
   human-branch substitution count (align each human HAR to its chimp ortholog
   — both coordinate sets are in the Cui supplement) is flagged as future work.

**Antagonistic pleiotropy.** HARs were deeply conserved across mammals, then
changed rapidly on the human lineage. When such an element sits in the
regulatory neighbourhood of a gene whose disruption causes a developmental
disorder *and* whose common variants associate with a psychiatric/cognitive
trait, the same regulatory tuning that may have been advantageous for human
brain evolution can raise disease liability under a different genetic
background — a hypothesis-generating pattern, not a mechanistic claim.

---

## License

MIT — see `LICENSE`.

## Citation of primary HAR source

Cui et al. (2025). *Nature* 640:991–999. doi:10.1038/s41586-025-08622-x. HAR set and neuronal
PLAC-seq HAR→gene interaction map.
