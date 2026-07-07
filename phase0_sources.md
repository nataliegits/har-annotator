# Phase 0 — Data source de-risking report

Genome build: **hg38** throughout. All sources below are already hg38-native — no
liftOver was required, so no intervals were dropped. (The pipeline retains a
liftOver hook for hg19 sources; it was not exercised.)

Every file is cached under `data/` and logged in `data/manifest.csv` with URL,
UTC access date, SHA256, and byte size.

## Core sources — all obtained

| # | Source | File | Rows | Notes |
|---|--------|------|------|-------|
| 1 | **HARs** — Cui et al. 2025, *Nature* (doi:10.1038/s41586-025-08622-x); coordinates redistributed in the `athenamarou/HAR-TFBS-Project` repo | `hars_hg38.bed` + `cui2025_HAR_supp4.xlsx` | 3,257 HARs | hg38 native. Union HAR set used in current neurodevelopmental HAR work; cross-referenced to HARsv2 / ZOOHAR IDs. |
| 2 | **Mammalian constraint** — Zoonomia 241-way phyloP | `cactus241way.phyloP.bw` (queried **remotely**, not downloaded) | per-base | UCSC `hg38/cactus241way`. bigWig is 9.6 GB; per-element remote query benchmarked at **~2 min for all 3,257 HARs** — the Phase-0 risk is cleared, so we use per-element phyloP mean/max directly (no fallback needed). |
| 3 | **Neurodevelopmental genes** — Gene2Phenotype Developmental Disorders (DDG2P) | `DDG2P_2026-06-28.csv.gz` | 2,860 gene–disease records | Gene symbol, confidence category, disease, allelic requirement, MONDO. |
| 4 | **Disease loci** — EBI GWAS Catalog, ontology-annotated associations | `gwas-catalog-associations_ontology-annotated-full.zip` | 1,150,105 associations | hg38 (`CHR_ID`/`CHR_POS`). ~33,400 rows match neuropsychiatric/neurodevelopmental traits (schizophrenia, autism, bipolar, ADHD, educational attainment, cognition, depression, epilepsy, intellectual disability). |
| 5 | **Developing-brain regulatory context** — ENCODE DNase-seq peaks, human brain embryo (105 days) | `ENCFF660HML_…_hg38.bed.gz` | 165,568 peaks | True mid-gestation cortex (embryonic life-stage), GRCh38. Used for "is this HAR an active regulatory element in the developing brain?" |

## Stretch sources — obtained

- **Enhancer→gene links (PLAC-seq in neurons):** the Cui 2025 supplement (sheet
  *"HARs interacting genes"*, 1,720 rows) provides HAR→interacting-gene links from
  neuronal PLAC-seq plus neuronal ATAC peak locations. This upgrades gene assignment
  beyond nearest-TSS for HARs that appear in the map.
- **Neuronal ATAC peaks:** embedded in the same supplement (promoter + distal
  interacting ATAC peak coordinates).

## Documented gaps / substitutions

- **SFARI Gene (autism) — UNAVAILABLE.** `gene-archive.sfari.org` is not
  reachable from this sandbox even after a network-access grant (the host is
  Cloudflare/exfil-gated at the proxy). The neurodevelopmental gene set is
  therefore **DDG2P alone** rather than SFARI ∪ DDG2P. DDG2P is a strong,
  clinically-curated developmental-disorder panel with confidence tiers, so the
  filter remains well-grounded; autism-specific genes are still represented via
  DDG2P entries. The pipeline is written so a SFARI CSV can be dropped in and
  unioned if it becomes available.
- **Acceleration score / human-substitution count — not in the coordinate
  source.** The Cui 2025 HAR table provides coordinates and chimp orthologs but
  not a per-HAR substitution count or acceleration score. Acceleration evidence
  in the scoring is therefore derived from **HAR width** (a proxy for
  substitution burden used cautiously) and the phyloP constraint signal, and this
  limitation is stated explicitly in the evidence schema. A substitution count
  can be added later by aligning human vs. chimp ortholog sequences (both
  coordinate sets are present).
- **Constraint fallback not needed.** The brief flagged the phyloP bigWig as the
  main risk; remote querying is fast, so per-element phyloP is used directly. The
  470-way phastCons constrained-element bigBed remains available at UCSC as a
  binary-overlap fallback if ever required.
