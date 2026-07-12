# HAR annotator — v2.0 (discovery mode)

![CI](https://github.com/nataliegits/har-annotator/actions/workflows/ci.yml/badge.svg)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)

> **This is the v2 companion to the shipped pipeline.** The main
> [`README.md`](../README.md) documents **v1**, the disease-anchored pipeline that
> requires each HAR to sit near a known neurodevelopmental *disease* gene. **v2
> drops that requirement.** Everything else — the constraint gate, the GWAS
> overlap, the seven scoring axes, the code — is identical. v1 is kept intact for
> reference; nothing here overwrites it.

**Human Accelerated Regions (HARs)** are stretches of the genome that stayed
frozen across mammalian evolution and then changed sharply on the human branch,
which makes them prime suspects for what rewired the human brain. v1 asks a
deliberately conservative question — *which HARs sit near a gene we already know
causes a brain disorder?* — and that disease-gene gate is doing real work, but it
can only ever return biology we have already named. **v2 opens the aperture:** it
keeps the same evolutionary and genetic-signal evidence but removes the
requirement that the nearby gene be a known disease gene, so the pipeline can
surface human-accelerated elements near genes with **no prior disease
annotation** at all.

Everything is hg38. Every data pull is cached, hashed, and logged to
`data/manifest.csv`, so the analysis replays end-to-end from source.

```
python demo_live_v2.py --discovery      # the v2 shortlist, live in <1 second
```

---

## What changes from v1 (and what does not)

| Funnel step | Constraint-only | v2 (discovery) | v1 (disease-anchored) |
|-------------|-----------------|----------------|-----------------------|
| 0 — All HARs (Cui et al. 2025, hg38) | 3,257 | 3,257 | 3,257 |
| A — Mammalian constraint (mean 241-way phyloP > 1.0) | **2,757** | 2,757 | 2,757 |
| B — Gene proximity (within 1 Mb of a TSS) | *not applied* | near **any** gene → 2,739 | near a **disease** gene → 1,718 |
| C — Within 25 kb of a neuropsychiatric GWAS lead SNP | *not applied* | **577** | **363** |

**Each level removes one gate.** Read the table right-to-left: v1 (disease-anchored,
**363**) is the strictest — constraint *and* GWAS overlap *and* a nearby DDG2P
disease gene. v2 (discovery, **577**) drops only the disease-gene requirement in
step B, keeping constraint and GWAS. **Constraint-only (2,757)** goes one gate
further and drops the GWAS overlap too, keeping only the evolutionary-constraint
filter — it is every HAR that passes step A, scored and ranked. The three sets
nest exactly: **363 ⊂ 577 ⊂ 2,757**.

**One caveat on the constraint-only level.** Because the acceleration and motif
axes need human–chimp alignments that were only fetched for the anchored/discovery
sets, the 2,757-element constraint-only shortlist is scored on **five axes**
(constraint, gene, disease, brain, temporal) with the two missing weights
redistributed proportionally (gene 0.25, disease 0.25, constraint 0.205, brain
0.148, temporal 0.148). Its ranking is therefore not directly comparable
element-for-element with the 7-axis anchored/discovery scores — treat it as a
wider net, not a finer ruler.

The result is a **strict superset**: 577 = the original **363** disease-anchored
elements + **214 newly surfaced** elements, every one of them near a gene that is
*not* in the disease panel.

---

## The three relaxation levels

The pipeline exposes three nested candidate sets, each answering a different
question. All three are precomputed tables you can load and re-rank; only the
first two have a live demo mode.

| Level | n | Command / file | The question it answers |
|-------|--:|----------------|-------------------------|
| **Disease-anchored** (v1) | 363 | `python run_pipeline.py` → `har_shortlist_ranked.parquet` | Which HARs sit near a gene we *already know* causes a brain disorder? |
| **Discovery** (v2) | 577 | `python demo_live_v2.py --discovery` → `har_shortlist_discovery.parquet` | Which constrained, GWAS-overlapping HARs sit near *any* gene, disease-annotated or not? |
| **Constraint-only** | 2,757 | `har_shortlist_relaxed.parquet` *(precomputed table — no demo flag)* | Which HARs are simply under strong mammalian constraint, before any disease/GWAS filter? |

The constraint-only table is scored on five axes (see the caveat above) and is
best used as a wide net for exploration. Re-rank it the same way as the discovery
set — a plain weighted sum over columns that already exist:

```python
import pandas as pd

df = pd.read_parquet("har_shortlist_relaxed.parquet")
# constraint-only ships 5 axes; weights renormalized to sum to 1
w = {"constraint": .205, "gene": .25, "disease": .25, "brain": .148, "temporal": .148}
df["total_score"] = sum(w[c] * df[f"score_{c}"] for c in w)
df = df.sort_values("total_score", ascending=False).reset_index(drop=True)
print(df.head(15)[["gene", "har_id", "total_score"]])
```

Relaxing to constraint-only reshuffles the leaderboard substantially: *TCF20*
rises to #1 and seven elements that the GWAS gate had excluded — *SMARCA2,
ITPR1, VLDLR, HOXC13, TRAPPC4, ORC4*, and a second *ZNF462* element — enter the
top 15. That is the point of the level: it shows what the disease/GWAS gates were
filtering out.

---

## The 214 new elements

These are the payload of v2: constrained, GWAS-overlapping human-accelerated
elements that v1 could never return because their nearest gene has no disease
annotation. The top of that new list reads like a who's-who of cortical
development that the disease-gene panel simply hasn't caught up to:

| HAR | Nearest gene | Why it's interesting |
|-----|--------------|----------------------|
| HAR_53   | **DAB1**  | reelin adaptor; controls neuronal migration / cortical layering |
| HAR_144  | **NPAS3** | bHLH TF; balanced-translocation disruption causes schizophrenia |
| HAR_3164 | **TLE4**  | Groucho co-repressor; deep-layer cortical neuron identity |
| HAR_3071 | **PTPRD** | synaptic adhesion phosphatase; ASD / ADHD associations |
| HAR_2425 | **EFNA5** | ephrin; axon guidance and topographic mapping |
| HAR_1073 | **CDH8**  | cadherin; autism-associated cell-adhesion |
| HAR_1767 / HAR_1361 | **TSHZ2 / TSHZ3** | teashirt zinc-finger developmental TFs |

Because these have no disease-gene tier, the `gene` axis (below) scores them on
proximity and PLAC-seq contact alone rather than on a disease-confidence tier, so
they generally rank *below* the disease-anchored hits — but they are now *in the
ranking*, visible and scored on identical axes, instead of filtered out before
they were ever seen.

Full 577-element output: `har_shortlist_discovery.parquet` / `.csv`. The 214
discovery-only elements on their own: `discovery_new_elements.csv`.

---

## What this pipeline does

v2 works in the same five moves as v1; only the gene test in move 1 differs:

1. **Filter.** Start from all 3,257 HARs; keep those that are evolutionarily
   constrained, sit **near any gene** (≤1 Mb), and overlap a genome-wide
   neuropsychiatric GWAS signal. 3,257 → **577**.
2. **Gather evidence.** For each survivor, record how strong every line of
   evidence is: constraint, gene link and distance, GWAS significance, fetal-brain
   activity, mid-fetal expression timing, human–chimp divergence, and TF-motif
   disruption.
3. **Score.** Collapse the evidence into one transparent weighted sum per HAR,
   every axis's contribution kept visible, and rank best-to-worst.
4. **Interpret.** Read the top of the list — now including the 214 new elements —
   and check it makes biological sense.
5. **Visualize & package.** Draw the figures (`fig_discovery_vs_anchored.png`)
   and bundle everything into a rerunnable package.

**One honest precision.** v2 still does not *prove* function. Dropping the
disease-gene gate widens what the pipeline can *see*; it does not change the fact
that the output is a **shortlist of candidate regulatory regions**, i.e.
hypotheses about which human-accelerated elements are worth studying. What v2 buys
is reach into un-annotated biology; what it costs is the disease-gene prior that
made v1's hits pre-vetted. Use v1 when you want conservative, disease-grounded
candidates; use v2 when you want to find something new.

---

## The transparent score (identical to v1)

The total is a **weighted sum of seven normalized (0–1) component scores**, no
black box. All 577 elements are scored together on the same scale, so the shared
363 and the 214 new are directly comparable.

```
total_score = Σ  WEIGHTS[c] · score_c        (c in the seven components below)
```

| Component | Default weight | Evidence |
|-----------|:-:|----------|
| `gene`         | 0.22 | disease-tier × TSS-distance decay (+0.15 if PLAC-seq-linked); **for the 214 new elements, no disease tier → proximity/PLAC only** |
| `disease`      | 0.22 | −log10(GWAS p) × proximity to lead SNP |
| `constraint`   | 0.18 | mean 241-way phyloP, min-max scaled across the scored set |
| `brain`        | 0.13 | ENCODE embryonic-brain DNase peak overlap (+ PLAC-seq ATAC), *where* |
| `temporal`     | 0.13 | fraction of the gene's prenatal expression in the mid-fetal window, *when* |
| `acceleration` | 0.07 | real human–chimp substitution rate (panTro5 ortholog alignment), min-max scaled |
| `motif`        | 0.05 | JASPAR TF-motif sites gained + lost between the human and chimp alleles, min-max scaled |

The only behavioural difference from v1 is in the `gene` axis: a disease gene
carries a confidence tier (definitive 1.0 → limited 0.3), whereas a v2-only gene
has no tier and is scored on distance decay and PLAC-seq contact alone. This is
why the disease-anchored 363 keep almost exactly their v1 order (Spearman
ρ = 0.9998, ZSWIM6 still #1) and the 214 new elements slot in below them.

For the full axis-by-axis math, weighting rationale, and caveats, see the main
[`README.md`](../README.md) — every axis is defined identically there.

---

## Install & run

```bash
# same conda env as v1 (Python 3.11)
conda create -n hars -c bioconda -c conda-forge python=3.11 \
    pandas numpy scipy matplotlib seaborn requests \
    pyranges=0.1.4 pybigwig biopython pyarrow openpyxl
conda activate hars

python run_pipeline.py                 # v1: reproduces the disease-anchored 363
python demo_live_v2.py --discovery     # v2: the 577-element discovery shortlist
python demo_live_v2.py --discovery --no-color   # strip ANSI for logs / projectors
```

The widest, constraint-only level (2,757 elements, GWAS gate dropped) is shipped
as a precomputed table rather than a demo mode — see
[**The three relaxation levels**](#the-three-relaxation-levels) below.

`demo_live_v2.py` is a strict superset of the shipped `demo_live.py`: it keeps
the three v1 moves (`default`, `--gene`, `--weight`) and adds `--discovery` (this
mode) and `--neglect W` (reward understudied target genes). v1's `demo_live.py`
is unchanged and remains the disease-anchored demo.

---

## Change the weights (the discovery shortlist re-ranks too)

v2 changes **which HARs enter** the shortlist (it drops the disease-gene gate);
it does **not** change how they are scored. The seven-axis weighted sum is
byte-identical to v1, and the shipped `har_shortlist_discovery.parquet` keeps
all seven `score_<axis>` columns — so you can re-rank the 577-element discovery
set under any weights, exactly the way v1's
[**Change the weights yourself**](README.md#change-the-weights-yourself)
section describes for the anchored 363.

The re-ranking is a plain weighted sum over columns that already exist, so it's
a few lines — no re-run, no network:

```python
import pandas as pd
from har_annotator.score import COMPONENTS   # the 7 axis names

df = pd.read_parquet("har_shortlist_discovery.parquet")
w  = {"constraint": .18, "gene": .10, "disease": .10,
      "brain": .13, "temporal": .13, "acceleration": .06, "motif": .30}  # your call
df["total_score"] = sum(w[c] * df[f"score_{c}"] for c in COMPONENTS)
df = df.sort_values("total_score", ascending=False).reset_index(drop=True)
df.to_csv("discovery_motif_weighted.csv", index=False)
```

Two things worth knowing, both the same as v1:

- **The `--weights` CLI flag lives on the *anchored* pipeline** (`python
  run_pipeline.py --weights motif=0.30`). It re-scores the disease-anchored 363.
  The discovery shortlist is a precomputed product; re-rank it with the snippet
  above (or edit the `WEIGHTS` dict in
  [`har_annotator/score.py`](har_annotator/score.py) and rebuild).
- **In the live demo the two moves are separate commands, not one.**
  `demo_live_v2.py --discovery` shows the discovery list; `demo_live_v2.py
  --weight motif=0.30` re-ranks the anchored 363 live. Discovery selects the
  candidate pool; `--weight`/the snippet sets how any pool is ranked — the two
  ideas are independent, they're just invoked separately in the shipped CLI.

---

## Interpretation & caveats

All of v1's caveats apply unchanged (overlap ≠ causation; nearest gene ≠ target
gene; acceleration is human–chimp divergence, not lineage-polarized; motif is a
model-based call). v2 adds one of its own:

- **No disease gate means no disease prior.** The 214 new elements have not been
  vetted against a curated disease-gene panel. Their gene links rest on proximity
  (and, where available, PLAC-seq contact) alone, so they are best read as
  *leads for follow-up*, not as pre-validated candidates. The disease-anchored 363
  remain the more conservative shortlist.

See `validation_ladder.md` for the wet-lab path (MPRA → Capture-C + CRISPRi →
humanized cortical organoids) that would test any of these candidates, v1 or v2.

---

## License

MIT (same as v1). Primary HAR source: Cui et al. 2025, *Nature* 640:991–999,
doi:10.1038/s41586-025-08622-x.
