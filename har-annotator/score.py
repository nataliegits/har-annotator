"""Phase 3 — transparent, inspectable evidence score.

The score is a simple weighted sum of six normalized component scores, each in
[0, 1]. Every component is retained as its own ``score_<component>`` column so
the total is fully decomposable — there is no black box. Weights are explicit
module constants and can be overridden per run.

    total_score = sum(WEIGHTS[c] * score_<c> for c in COMPONENTS)

Components
----------
constraint     : mammalian constraint strength (phyloP mean), min-max scaled
acceleration   : substitution-burden proxy (HAR width), log-scaled  [see phase0 caveat]
gene           : neurodev-gene link quality — DDG2P confidence tier, boosted if
                 the gene is assigned by Hi-C rather than nearest-TSS, and
                 decaying with HAR–TSS distance
disease        : neuropsychiatric GWAS overlap — significance (−log10 p) and
                 proximity to the lead SNP
brain          : active in the developing brain — ENCODE embryonic-cortex DNase
                 peak overlap (+ Hi-C ATAC support)
motif          : reserved (0 unless a motif-disruption annotation is supplied)
"""
from __future__ import annotations

import numpy as np
import pandas as pd

COMPONENTS = ["constraint", "acceleration", "gene", "disease", "brain", "motif"]

# Default weights — deliberately simple and legible. Constraint + disease +
# gene carry the most weight (they define the biological question); brain
# activity and acceleration are supporting; motif is reserved.
WEIGHTS = {
    "constraint": 0.20,
    "acceleration": 0.10,
    "gene": 0.25,
    "disease": 0.25,
    "brain": 0.15,
    "motif": 0.05,
}

_CONF_SCORE = {"definitive": 1.0, "strong": 0.75, "moderate": 0.5,
               "limited": 0.3, "disputed": 0.1, "refuted": 0.0}


def _minmax(x: pd.Series) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce")
    lo, hi = x.min(), x.max()
    if not np.isfinite(lo) or hi == lo:
        return pd.Series(np.zeros(len(x)), index=x.index)
    return (x - lo) / (hi - lo)


def compute_scores(ev: pd.DataFrame, weights: dict | None = None,
                   gene_window: int = 1_000_000,
                   gwas_window: int = 25_000) -> pd.DataFrame:
    """Return ``ev`` with per-component score columns and a ``total_score``,
    sorted descending. All component columns are in [0, 1]."""
    w = {**WEIGHTS, **(weights or {})}
    df = ev.copy()

    # constraint: min-max scaled mean phyloP
    df["score_constraint"] = _minmax(df["phylop_mean"])

    # acceleration proxy: log HAR width, min-max scaled
    df["score_acceleration"] = _minmax(np.log1p(df["width"]))

    # gene: confidence tier * distance-decay, +0.15 boost for Hi-C assignment
    conf = df["gene_confidence"].str.lower().map(_CONF_SCORE).fillna(0.3)
    dist_decay = 1.0 - (df["gene_distance"].clip(0, gene_window) / gene_window)
    hic_boost = np.where(df["gene_assignment_method"].eq("hic_linked"), 0.15, 0.0)
    df["score_gene"] = np.clip(conf * dist_decay + hic_boost, 0, 1)

    # disease: significance (−log10 p, capped at 50) * proximity decay
    neglogp = -np.log10(pd.to_numeric(df["gwas_pval"], errors="coerce").clip(lower=1e-300))
    sig = (neglogp.clip(0, 50) / 50.0)
    prox = 1.0 - (df["gwas_distance"].clip(0, gwas_window) / gwas_window)
    df["score_disease"] = np.clip(sig.fillna(0) * prox.fillna(0), 0, 1)

    # brain: DNase peak overlap (0.8) + Hi-C ATAC support (0.2)
    df["score_brain"] = (df["brain_dnase_overlap"].astype(float) * 0.8
                         + df["hic_atac"].fillna(False).astype(float) * 0.2)

    # motif: reserved
    if "score_motif" not in df:
        df["score_motif"] = df.get("motif_disruption", pd.Series(0.0, index=df.index)).fillna(0.0)

    contrib = []
    for c in COMPONENTS:
        col = f"contrib_{c}"
        df[col] = w[c] * df[f"score_{c}"]
        contrib.append(col)
    df["total_score"] = df[contrib].sum(axis=1)

    df = df.sort_values("total_score", ascending=False).reset_index(drop=True)
    df.insert(0, "rank", np.arange(1, len(df) + 1))
    return df
