"""
Compute allelic polarization results for human-chimp (HC) and human-gorilla (HG)
hybrid expression experiments.

Outputs
-------
HC: outputs_05Jan2026_humanMPRA_draft1/ExpLBM_polarization_results.tsv
HG: outputs_24May2026_humanMPRA_draft2/ExpLBM_HG_polarization_results.tsv

Usage
-----
python polarization.py          # run both pipelines
python polarization.py --hc     # run only the HC pipeline
python polarization.py --hg     # run only the HG pipeline
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


# ── PATHS ──────────────────────────────────────────────────────────────────────

_HC_BASE = Path(
    "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids"
    "/human_chimp/ExpLBM/outputs_05Jan2026_humanMPRA_draft1"
)
_HG_BASE = Path(
    "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids"
    "/human_gorilla/ExpLBM/outputs_05Jan2026_humanMPRA_draft1"
)

HC_ASE_INFO = Path(
    "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids"
    "/human_chimp/summarized_data/ASE_info.tsv"
)

HC_DESEQ2_HUMAN_REF  = _HC_BASE / "ASE/DESeq2_aligned2human.txt"
HC_DESEQ2_CHIMP_REF  = _HC_BASE / "ASE/DESeq2_aligned2chimp.txt"

HG_DESEQ2_HUMAN_REF  = _HG_BASE / "ASE/DESeq2_aligned2human.txt"

HG_ASE_METADATA      = _HG_BASE / "ASE/combined_lfc_with_ASE_definition.tsv"

HC_OUTPUT = _HC_BASE / "ExpLBM_polarization_results.tsv"
HG_OUTPUT = Path(
    "/home/labs/davidgo/Collaboration/USEFUL_DATASETS/Expression/Hybrids"
    "/human_gorilla/ExpLBM/outputs_24May2026_humanMPRA_draft2"
    "/ExpLBM_HG_polarization_results.tsv"
)


# ── SHARED UTILITIES ───────────────────────────────────────────────────────────

def _log2tpm1(series: pd.Series) -> pd.Series:
    return np.log2(series.astype(float) + 1)


# ── HC PIPELINE ────────────────────────────────────────────────────────────────

def load_hc_data() -> pd.DataFrame:
    """
    Load HC ASE info and gorilla TPMs from HG experiment.

    Returns a DataFrame indexed by gene with:
    - ExpLBM_* columns from ASE_info.tsv (human/chimp TPMs, LFCs, ASE type)
    - mean_gorilla_tpm_HG, mean_human_tpm_HG, mean_total_tpm_HG from HG DESeq2
    """
    ase_table = (
        pd.read_csv(HC_ASE_INFO, sep="\t", dtype={9: str})
        .set_index("Gene")
        .filter(regex="^ExpLBM", axis=1)
    )

    hg_tpm = (
        pd.read_csv(HG_DESEQ2_HUMAN_REF, sep="\t", index_col=0)
        [["mean_ref_tpm", "mean_alt_tpm", "mean_total_tpm"]]
        .rename(columns={
            "mean_ref_tpm":   "mean_human_tpm_HG",
            "mean_alt_tpm":   "mean_gorilla_tpm_HG",
            "mean_total_tpm": "mean_total_tpm_HG",
        })
    )

    return ase_table.join(hg_tpm, how="left")


def compute_hc_polarization(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute polarization score and derived classification for HC experiment.

    Polarization logic
    ------------------
    - G_clipped: gorilla TPM clipped to [min(H, C), max(H, C)]
    - human_derived_score = 2 * ((H - G_clipped) / (H - C) - 0.5)  in [-1, 1]
    - Distances in log2 space: d_GH = |log2G - log2H|, d_GC = |log2G - log2C|
    - derived label (ASE genes only):
        - "human-derived" if d_GC < d_GH  (gorilla closer to chimp)
        - "chimp-derived" if d_GH < d_GC  (gorilla closer to human)
    """
    df = df.copy()

    H = df["ExpLBM_TPM_human_allele"].astype(float)
    C = df["ExpLBM_TPM_chimp_allele"].astype(float)
    G = df["mean_gorilla_tpm_HG"].astype(float)

    G_clip = G.clip(lower=np.minimum(H, C), upper=np.maximum(H, C))
    den = H - C
    score = 2 * (((H - G_clip) / den) - 0.5)
    df["G_clipped"] = G_clip
    df["human_derived_score"] = score.where(den != 0, np.nan)

    df["ExpLBM_log2TPM_gorilla_allele"]    = _log2tpm1(df["mean_gorilla_tpm_HG"])
    df["ExpLBM_log2TPM_human_allele"]      = _log2tpm1(df["ExpLBM_TPM_human_allele"])
    df["ExpLBM_log2TPM_chimp_allele"]      = _log2tpm1(df["ExpLBM_TPM_chimp_allele"])
    df["ExpLBM_log2TPM_human_allele_HG"]   = _log2tpm1(df["mean_human_tpm_HG"])
    df["ExpLBM_log2TPM_total"]             = _log2tpm1(df["ExpLBM_TPM_total"])

    df.rename(columns={
        "mean_gorilla_tpm_HG": "ExpLBM_TPM_gorilla_allele",
        "mean_human_tpm_HG":   "ExpLBM_TPM_human_allele_HG",
    }, inplace=True)

    df["d_GH"] = np.abs(df["ExpLBM_log2TPM_gorilla_allele"] - df["ExpLBM_log2TPM_human_allele"])
    df["d_GC"] = np.abs(df["ExpLBM_log2TPM_gorilla_allele"] - df["ExpLBM_log2TPM_chimp_allele"])

    is_ASE    = df["ExpLBM_gene_ase_type"].notna()
    has_gorilla = df["ExpLBM_log2TPM_gorilla_allele"].notna()

    df["derived"] = pd.Series(pd.NA, index=df.index, dtype="object")
    df.loc[is_ASE & has_gorilla & (df["d_GH"] < df["d_GC"]), "derived"] = "chimp-derived"
    df.loc[is_ASE & has_gorilla & (df["d_GC"] < df["d_GH"]), "derived"] = "human-derived"

    return df[[
        "ExpLBM_TPM_human_allele",
        "ExpLBM_TPM_human_allele_HG",
        "ExpLBM_TPM_chimp_allele",
        "ExpLBM_TPM_gorilla_allele",
        "ExpLBM_TPM_total",
        "ExpLBM_log2TPM_human_allele",
        "ExpLBM_log2TPM_human_allele_HG",
        "ExpLBM_log2TPM_chimp_allele",
        "ExpLBM_log2TPM_gorilla_allele",
        "ExpLBM_log2TPM_total",
        "ExpLBM_LFC_human_ref",
        "ExpLBM_LFC_pvalue_human_ref",
        "ExpLBM_LFC_padj_human_ref",
        "ExpLBM_LFC_chimp_ref",
        "ExpLBM_LFC_pvalue_chimp_ref",
        "ExpLBM_LFC_padj_chimp_ref",
        "ExpLBM_gene_ase_type",
        "d_GH",
        "d_GC",
        "G_clipped",
        "human_derived_score",
        "derived",
    ]]


def run_hc_pipeline() -> None:
    print("Running HC polarization pipeline...")
    df = load_hc_data()
    result = compute_hc_polarization(df)
    HC_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(HC_OUTPUT, sep="\t")
    print(f"  Saved {len(result)} genes → {HC_OUTPUT}")


# ── HG PIPELINE ────────────────────────────────────────────────────────────────

def load_hg_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Load all inputs for the HG polarization pipeline.

    Returns
    -------
    hc_chimp_tpm : DataFrame
        HC chimp allele TPM (from ASE_info.tsv), indexed by gene.
    hg_tpm : DataFrame
        HG human/gorilla allele TPMs and totals (from HG DESeq2 aligned to human).
    ase_hg : DataFrame
        HG ASE metadata: LFCs, p-values, ASE definition (from combined_lfc table).
    """
    hc_chimp_tpm = (
        pd.read_csv(HC_ASE_INFO, sep="\t", dtype={9: str})
        .set_index("Gene")
        [["ExpLBM_TPM_chimp_allele"]]
    )

    hg_tpm = (
        pd.read_csv(HG_DESEQ2_HUMAN_REF, sep="\t", index_col=0)
        [["mean_ref_tpm", "mean_alt_tpm", "mean_total_tpm"]]
        .rename(columns={
            "mean_ref_tpm":   "ExpLBM_TPM_human_allele_HG",
            "mean_alt_tpm":   "ExpLBM_TPM_gorilla_allele_HG",
            "mean_total_tpm": "ExpLBM_TPM_total_HG",
        })
    )

    ase_hg = pd.read_csv(HG_ASE_METADATA, sep="\t", index_col=0)

    return hc_chimp_tpm, hg_tpm, ase_hg


def compute_hg_polarization(
    hc_chimp_tpm: pd.DataFrame,
    hg_tpm: pd.DataFrame,
    ase_hg: pd.DataFrame,
) -> pd.DataFrame:
    """
    Compute polarization distances and classification for HG experiment.

    Polarization logic
    ------------------
    Uses chimp (from HC) as outgroup reference for human vs gorilla comparison:
    - Distances in log2 space: d_CH = |log2C_HC - log2H_HG|, d_CG = |log2C_HC - log2G_HG|
    - derived label (all genes with chimp data):
        - "human-derived"  if d_CG < d_CH  (chimp closer to gorilla)
        - "gorilla-derived" if d_CH < d_CG  (chimp closer to human)
    """
    df = (
        hc_chimp_tpm
        .rename(columns={"ExpLBM_TPM_chimp_allele": "ExpLBM_TPM_chimp_allele_HC"})
        .join(hg_tpm, how="left")  # keep all ASE genes; genes with no HG data get NaN
    )

    df["ExpLBM_log2TPM_chimp_allele_HC"]   = _log2tpm1(df["ExpLBM_TPM_chimp_allele_HC"])
    df["ExpLBM_log2TPM_human_allele_HG"]   = _log2tpm1(df["ExpLBM_TPM_human_allele_HG"])
    df["ExpLBM_log2TPM_gorilla_allele_HG"] = _log2tpm1(df["ExpLBM_TPM_gorilla_allele_HG"])
    df["ExpLBM_log2TPM_total_HG"]          = _log2tpm1(df["ExpLBM_TPM_total_HG"])

    df["d_CH"] = np.abs(df["ExpLBM_log2TPM_chimp_allele_HC"] - df["ExpLBM_log2TPM_human_allele_HG"])
    df["d_CG"] = np.abs(df["ExpLBM_log2TPM_chimp_allele_HC"] - df["ExpLBM_log2TPM_gorilla_allele_HG"])

    has_chimp = df["ExpLBM_log2TPM_chimp_allele_HC"].notna()
    df["derived"] = pd.Series(pd.NA, index=df.index, dtype="object")
    df.loc[has_chimp & (df["d_CH"] < df["d_CG"]), "derived"] = "gorilla-derived"
    df.loc[has_chimp & (df["d_CG"] < df["d_CH"]), "derived"] = "human-derived"

    lfc_cols = ase_hg[["log2FoldChange_human_ref", "pvalue_human_ref", "padj_human_ref", "ASE_definition"]].rename(
        columns={
            "log2FoldChange_human_ref": "HG_LFC_human_ref",
            "pvalue_human_ref":         "HG_pvalue_human_ref",
            "padj_human_ref":           "HG_padj_human_ref",
        }
    )
    df = df.join(lfc_cols, how="left")

    return df[[
        "ExpLBM_TPM_human_allele_HG",
        "ExpLBM_TPM_chimp_allele_HC",
        "ExpLBM_TPM_gorilla_allele_HG",
        "ExpLBM_TPM_total_HG",
        "ExpLBM_log2TPM_human_allele_HG",
        "ExpLBM_log2TPM_chimp_allele_HC",
        "ExpLBM_log2TPM_gorilla_allele_HG",
        "ExpLBM_log2TPM_total_HG",
        "HG_LFC_human_ref",
        "HG_pvalue_human_ref",
        "HG_padj_human_ref",
        "ASE_definition",
        "derived",
        "d_CH",
        "d_CG",
    ]]


def run_hg_pipeline() -> None:
    print("Running HG polarization pipeline...")
    hc_chimp_tpm, hg_tpm, ase_hg = load_hg_data()
    result = compute_hg_polarization(hc_chimp_tpm, hg_tpm, ase_hg)
    HG_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(HG_OUTPUT, sep="\t")
    print(f"  Saved {len(result)} genes → {HG_OUTPUT}")


# ── MAIN ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute allelic polarization results for HC and/or HG hybrids."
    )
    parser.add_argument("--hc", action="store_true", help="Run human-chimp pipeline only")
    parser.add_argument("--hg", action="store_true", help="Run human-gorilla pipeline only")
    args = parser.parse_args()

    run_both = not (args.hc or args.hg)

    if args.hc or run_both:
        run_hc_pipeline()
    if args.hg or run_both:
        run_hg_pipeline()


if __name__ == "__main__":
    main()
