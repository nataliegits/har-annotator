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
**HAR_575 → *POC1B*** (Hi-C-linked, educational-attainment GWAS overlap,
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
| `gene`         | 0.25 | DDG2P confidence tier × TSS-distance decay, +0.15 if Hi-C-linked |
| `disease`      | 0.25 | −log10(GWAS p) × proximity to lead SNP |
| `brain`        | 0.15 | ENCODE embryonic-brain DNase peak overlap (+ Hi-C ATAC support) |
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

| Source | Used for | Access |
|--------|----------|--------|
| **HARs** — Cui et al. 2025, *Nature* (doi:10.1038/s41586-025-08622), Supp. Table 4 | 3,257 HARs, hg38-native; neuronal Hi-C HAR→gene links | GitHub mirror `athenamarou/HAR-TFBS-Project` |
| **Constraint** — Zoonomia 241-way phyloP bigWig (UCSC `cactus241way`) | per-HAR mean/max phyloP | **queried remotely** (9.6 GB; not downloaded) |
| **Neurodev genes** — DDG2P (Gene2Phenotype Developmental Disorders) | 2,524 genes mapped to hg38 coords, with confidence tier | EBI FTP |
| **GWAS** — EBI GWAS Catalog (ontology-annotated) | 22,489 genome-wide-sig neuropsychiatric lead SNPs | EBI FTP |
| **Developing brain** — ENCODE embryonic-brain DNase-seq (`ENCFF660HML`, male embryo 105 days, expt ENCSR420RWU) | 165,568 peaks; regulatory-activity axis | ENCODE portal |
| **Gene coords** — UCSC refGene hg38 | gene symbol → strand-aware TSS | UCSC |

**Substitution note.** SFARI Gene was unreachable from the build environment
(`gene-archive.sfari.org` blocked), so the neurodevelopmental gene set is
**DDG2P only**. The pipeline still unions any SFARI symbols supplied via
`--sfari`, so a SFARI export can be dropped in without code changes.

---

## Package layout

```
har_annotator/
  download.py    cached + hashed + manifest-logged fetches
  data_io.py     load HAR base table + Cui Hi-C sheet from the supplement
  references.py  build derived tables: genes, neurodev (DDG2P), GWAS loci, Hi-C links
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

## Interpretation & caveats

Read `top_candidates_interpretation.md` for the top-candidate biology. Three
caveats apply throughout and are load-bearing:

1. **Overlap ≠ causation.** "Overlaps a GWAS locus" means within 25 kb of a
   lead SNP — not that the HAR is the causal element or that its
   human-specific changes drive the association.
2. **Nearest gene ≠ target gene.** `nearest_tss` assignments are a proximity
   heuristic; `hic_linked` assignments rest on a neuronal Hi-C loop (stronger,
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

Cui et al. (2025). *Nature*. doi:10.1038/s41586-025-08622. HAR set and neuronal
Hi-C HAR→gene interaction map.
