import argparse
import sys
import pandas as pd
import numpy as np
from scipy.stats import ranksums
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import os

sys.path.insert(0, '/home/labs/davidgo/Collaboration/backup/MPRA_QC_pipeline/QC_analysis/scripts')
import const

# python ase_sliding_window.py \
#     --gff genes.gff3 \
#     --counts counts.tsv \
#     --window 50 \
#     --out ase_output.tsv

# Standard font size and family
FONT_SIZE_small = 18 #changed from 16 NM 06-02-2025
FONT_SIZE_big = 20 #changed from 18 NM 06-02-2025

FONT_FAMILY = 'sans-serif'

plt.rcParams.update({
    #'axes.titlesize': FONT_SIZE_big,
    #'axes.labelsize': FONT_SIZE_big,
    #'xtick.labelsize': FONT_SIZE_small,
    #'ytick.labelsize': FONT_SIZE_small,
    #'legend.fontsize': FONT_SIZE_small,
    #'legend.title_fontsize': FONT_SIZE_small,
    'axes.labelweight': 'bold',
    'font.family': FONT_FAMILY, 
    'font.sans-serif': ['Arial'], 
    #'axes.linewidth': 1.0,
    #'figure.figsize': (8, 8),
    #'axes.grid': False,  # No grid
    'axes.spines.top': False,  # Top border off
    'axes.spines.right': False,  # Right border off
    #'figure.facecolor': 'none',
    #'axes.facecolor': 'none',  
#        'legend.frameon': False  # No frame for legend
})



# --------------------------------------------------------------------
# Parse arguments
# --------------------------------------------------------------------
def get_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True, help="GFF file with gene coordinates")
    ap.add_argument("--counts", required=True, help="TSV with gene, ref_counts, alt_counts")
    ap.add_argument("--window", type=int, default=50, help="Sliding window size in genes")
    ap.add_argument("--out", required=True, help="Output TSV file")
    return ap.parse_args()

def load_gtf(gtf_file):
    """
    Minimal GTF parser:
    - keeps only rows with feature == 'gene'
    - extracts gene_name (preferred) or gene_id
    - returns DataFrame with columns: gene, chrom, position
    """
    rows = []
    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            chrom, source, feature, start, end, score, strand, phase, attributes = parts

            # Only gene features
            if feature != "gene":
                continue

            # Parse GTF-style attributes: key "value";
            attr_dict = {}
            for field in attributes.split(";"):
                field = field.strip()
                if not field:
                    continue
                # Expect: key "value"
                if " " not in field:
                    continue
                key, val = field.split(" ", 1)
                key = key.strip()
                val = val.strip().strip('"')  # remove quotes
                attr_dict[key] = val

            gene_name = attr_dict.get("gene_name") or attr_dict.get("gene_id")

            rows.append([gene_name, chrom, int(start)])

    df = pd.DataFrame(rows, columns=["gene", "chrom", "position"])
    return df

# --------------------------------------------------------------------
# Compute ASE
# --------------------------------------------------------------------
def compute_ase(df, min_counts=10):
    # total counts
    df["total"] = df["ref_counts"] + df["alt_counts"]

    # prepare columns for safe division
    df["ref2"] = df["ref_counts"].replace(0, 1)
    df["alt2"] = df["alt_counts"].replace(0, 1)

    # compute ase
    df["ase"] = np.log2(df["ref2"] / df["alt2"])

    # apply low-count filter: ASE = NaN when total <= 9
    df.loc[df["total"] < min_counts, "ase"] = np.nan

    # cleanup
    return df.drop(columns=["ref2", "alt2", "total"])

# --------------------------------------------------------------------
# Sliding window median and Wilcoxon
# --------------------------------------------------------------------
def sliding_stats(ase_vals, global_ase, window):
    """
    ase_vals   : ASE values for this chromosome (1D array-like)
    global_ase : ASE values for the whole genome (df['ase'].values)
    window     : window size in genes

    Returns:
      med   : median ASE per window (length n_win)
      pvals : Wilcoxon p-value per window (length n_win)
    """
    ase_vals = np.asarray(ase_vals)
    global_ase = np.asarray(global_ase)

    n = len(ase_vals)
    if n < window:
        return np.array([]), np.array([])

    n_win = n - window + 1
    med   = np.full(n_win, np.nan)
    pvals = np.full(n_win, np.nan)

    from scipy.stats import ranksums

    # Background = all finite ASE values genome-wide
    bg_all = global_ase[np.isfinite(global_ase)]
    if len(bg_all) == 0:
        return med, pvals  # all NaN

    for i in range(n_win):
        # if the gene at the start of the window has no ASE, skip
        if not np.isfinite(ase_vals[i]):
            continue

        win = ase_vals[i:i+window]
        win_ok = win[np.isfinite(win)]
        if len(win_ok) == 0:
            continue

        med[i] = np.median(win_ok)
        try:
            pvals[i] = ranksums(win_ok, bg_all, alternative="two-sided").pvalue
        except Exception:
            pvals[i] = np.nan

    return med, pvals

def explain_row(df, chrom, gene, window):
    sub = df[df["chrom"] == chrom].reset_index(drop=True)
    i = sub.index[sub["gene"] == gene][0]

    win = sub.loc[i : i + window - 1, ["gene", "ase"]]
    ase_vals = win["ase"].values
    ase_ok   = ase_vals[np.isfinite(ase_vals)]
    bg = np.concatenate([
        sub.loc[: i - 1, "ase"].values,
        sub.loc[i + window :, "ase"].values
    ])
    bg_ok = bg[np.isfinite(bg)]

    p = ranksums(ase_ok, bg_ok).pvalue

    print(f"Chrom {chrom}, gene {gene}, index {i}")
    print("\nWindow genes & ASE:")
    print(win)
    print("\nMedian (recalc vs stored):", np.median(ase_ok), sub.loc[i, "median_window"])
    print("Wilcoxon p (recalc vs stored):", p, sub.loc[i, "wilcoxon_p"])


def chrom_sort_key(c):
    c = str(c)
    if c.startswith("chr"):
        c = c[3:]  # strip "chr"
    # numeric chromosomes
    if c.isdigit():
        return (0, int(c))
    # special ones
    mapping = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return (1, mapping.get(c, 999))  # unknown stuff at the end


def plot_medians_by_chrom(df, window, out_file="median_by_chromosomes.pdf"):

    # Extract chromosomes in natural order (chr1, chr2, …)
    chroms = sorted(df["chrom"].unique(), key=chrom_sort_key)

    n_chr = len(chroms)
    n_row = int(np.ceil(np.sqrt(n_chr)))
    n_col = int(np.ceil(n_chr / n_row))

    # const.set_plot_style()
    fig, axes = plt.subplots(n_row, n_col, figsize=(14, 12), sharey=False)
    axes = axes.flatten()

    for i, (ax, chrom) in enumerate(zip(axes, chroms)):
        print(chr)

        sub = df[df["chrom"] == chrom].reset_index(drop=True)

        # Prepare x-axis positions for the medians
        # Use the genomic position of the gene at the start of each window
        n = len(sub)
        if n < window:
            ax.set_title(f"{chrom} (too few genes)", fontsize=16)
            ax.axis("off")
            continue

        x_pos = sub.loc[: n - window, "position"]           # gene positions for windows
        y_med = sub.loc[: n - window, "median_window"]      # aligned medians

        # Remove NaNs (from low-coverage)
        mask = np.isfinite(y_med)
        x_pos = x_pos[mask]
        y_med = y_med[mask]

        ax.scatter(x_pos, y_med, s=8)
        ax.axhline(0, color="gray", linestyle="--", linewidth=0.7)

        ax.set_title(chrom, fontsize=16)

        x = np.asarray(x_pos, dtype=float)
        x = x[np.isfinite(x)]

        if x.size:
            xmin, xmax = x.min(), x.max()
            ax.set_xticks([xmin, xmax])
            ax.xaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
        else:
            ax.set_xticks([])

        # xticks = ax.get_xticks()
        # ax.set_xticks([xticks[0], xticks[-1]])
        ax.set_ylim(-4, 4)

        if i == 0:
            ax.set_xlabel("Relative Genomic Position", fontsize=14)
            ax.set_ylabel(r"Median $\log_{2}$(human/chimp)", fontsize=14)
            ax.yaxis.label.set_multialignment('center')  # nice centering for multi-line
            ax.tick_params(axis="y", left=True, labelleft=True)
        else:
            ax.set_xlabel("")
            ax.set_ylabel("")
            ax.tick_params(axis="y", left=False, labelleft=False)


    # Hide unused axes if any
    for ax in axes[len(chroms):]:
        ax.axis("off")

    plt.tight_layout()
    fig.savefig(f"{out_file}.png", dpi=500,bbox_inches='tight',transparent=True) 
    fig.savefig(f"{out_file}.svg", bbox_inches="tight", transparent=True)
    fig.savefig(f"{out_file}.eps", dpi=500,bbox_inches='tight') # increasd DPI to 500 NM 19/09
    fig.savefig(f"{out_file}.pdf", dpi=500,bbox_inches='tight',transparent=True) 
    #plt.savefig(f"{out_file}.pdf")
    plt.close()

    print(f"Saved median plot to {out_file}.pdf")

def plot_median_and_p_per_chrom(df, window, out_prefix="chr_medians_pvals"):
    """
    For each chromosome in df, create a PDF:
      top:  median_window vs genomic position
      bottom: -log10(wilcoxon_p) vs genomic position
    """

    chroms = sorted(df["chrom"].unique(), key=chrom_sort_key)

    # Global Bonferroni threshold like in your R code:
    # bonf = -log10(0.05 / length(ase))
    n_ase = df["ase"].notna().sum()
    if n_ase > 0:
        bonf = -np.log10(0.05 / n_ase)
    else:
        bonf = None

    for chrom in chroms:
        sub = df[df["chrom"] == chrom].reset_index(drop=True)

        # If chromosome has fewer genes than window, just skip/annotate
        if len(sub) < window:
            print(f"Skipping {chrom}: fewer than {window} genes")
            continue

        # --- median data ---
        med = sub["median_window"].values
        pos = sub["position"].values
        mask_med = np.isfinite(med)

        # --- p-values ---
        pvals = sub["wilcoxon_p"].values
        with np.errstate(divide="ignore", invalid="ignore"):
            stats_log = -np.log10(pvals)
        mask_p = np.isfinite(stats_log)

        # If no valid medians or p-values, skip
        if not mask_med.any() and not mask_p.any():
            print(f"Skipping {chrom}: no finite medians / p-values")
            continue

        fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
        ax1, ax2 = axes

        # ---------------------------
        # 1) Median window plot
        # ---------------------------
        if mask_med.any():
            ax1.scatter(pos[mask_med], med[mask_med], s=8)
        ax1.axhline(0, linestyle="--", linewidth=0.7)
        ax1.set_ylabel("Median log2(ASE)")
        ax1.set_title(f"{chrom} - window medians (window = {window} genes)")

        # ---------------------------
        # 2) Wilcoxon p-value plot
        # ---------------------------
        if mask_p.any():
            ax2.scatter(pos[mask_p], stats_log[mask_p], s=8)
        if bonf is not None:
            ax2.axhline(bonf, linestyle="--", linewidth=0.7)
        ax2.set_ylabel("-log10(p)")
        ax2.set_xlabel("Genomic position")

        plt.tight_layout()
        out_file = f"{out_prefix}_{chrom}.pdf"
        # save pvals
        pval_out_file = f"{out_prefix}_{chrom}_pvals.tsv"
        plt.savefig(out_file)
        plt.close()
        print(f"Saved {out_file}")


# ---------- Using TPM ----------

def compute_log2_tpm(df, tpm_col="ref_tpm", out_col="log2_ref_tpm", pseudocount=0.1):
    """
    Add log2(TPM + pseudocount) column.
    pseudocount helps avoid -inf for TPM=0.
    """
    if tpm_col not in df.columns:
        raise KeyError(f"Missing column '{tpm_col}' in df.")
    vals = pd.to_numeric(df[tpm_col], errors="coerce")
    df[out_col] = np.log2(vals + pseudocount)
    return df

def sliding_median(values, window):
    """
    Sliding median for a 1D array, returns length n-window+1.
    NaNs are ignored within each window (median of finite values).
    If a window has no finite values -> NaN.
    """
    values = np.asarray(values, dtype=float)
    n = len(values)
    if n < window:
        return np.array([], dtype=float)

    out = np.empty(n - window + 1, dtype=float)
    for i in range(n - window + 1):
        w = values[i:i+window]
        w = w[np.isfinite(w)]
        out[i] = np.median(w) if w.size else np.nan
    return out

def plot_sliding_medians_by_chrom_tpm(
    df,
    window,
    tpm_log_col="log2_ref_tpm",
    out_file="median_by_chromosomes_refTPM",
    y_label=r"Median $\log_{2}$(ref TPM + pc)",
    y_lim=None,
    point_size=8,
):
    """
    Like plot_medians_by_chrom, but uses log2(TPM) sliding medians per chromosome.
    Expects df sorted or at least has chrom + position columns.
    """

    if "chrom" not in df.columns or "position" not in df.columns:
        raise KeyError("df must contain 'chrom' and 'position' columns.")
    if tpm_log_col not in df.columns:
        raise KeyError(f"df must contain '{tpm_log_col}' (run compute_log2_tpm first).")

    # natural order helper: assumes you already have chrom_sort_key in your script
    chroms = sorted(df["chrom"].dropna().unique(), key=chrom_sort_key)

    n_chr = len(chroms)
    n_row = int(np.ceil(np.sqrt(n_chr)))
    n_col = int(np.ceil(n_chr / n_row))

    fig, axes = plt.subplots(n_row, n_col, figsize=(14, 12), sharey=False)
    axes = np.array(axes).flatten()

    for i, (ax, chrom) in enumerate(zip(axes, chroms)):
        sub = df[df["chrom"] == chrom].sort_values("position").reset_index(drop=True)
        n = len(sub)

        if n < window:
            ax.set_title(f"{chrom} (too few genes)", fontsize=16)
            ax.axis("off")
            continue

        vals = sub[tpm_log_col].to_numpy(dtype=float)
        med = sliding_median(vals, window)  # length n-window+1

        # x = position of the *start* of each window (same convention as your ASE plot)
        x_pos = sub.loc[: n - window, "position"].to_numpy(dtype=float)

        # filter NaNs
        mask = np.isfinite(med) & np.isfinite(x_pos)
        x_pos = x_pos[mask]
        med = med[mask]

        ax.scatter(x_pos, med, s=point_size)
        ax.set_title(chrom, fontsize=16)

        # x ticks: min/max
        if x_pos.size:
            xmin, xmax = float(np.min(x_pos)), float(np.max(x_pos))
            ax.set_xticks([xmin, xmax])
            ax.xaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
        else:
            ax.set_xticks([])

        if y_lim is not None:
            ax.set_ylim(y_lim)

        # label only first panel
        if i == 0:
            ax.set_xlabel("Genomic Position", fontsize=14)
            ax.set_ylabel(y_label, fontsize=14)
            ax.yaxis.label.set_multialignment('center')
            ax.tick_params(axis="y", left=True, labelleft=True)
        else:
            ax.set_xlabel("")
            ax.set_ylabel("")
            ax.tick_params(axis="y", left=False, labelleft=False)

    # Hide unused axes
    for ax in axes[len(chroms):]:
        ax.axis("off")

    plt.tight_layout()
    fig.savefig(f"{out_file}.png", dpi=500, bbox_inches="tight", transparent=True)
    fig.savefig(f"{out_file}.svg", bbox_inches="tight", transparent=True)
    fig.savefig(f"{out_file}.eps", dpi=500, bbox_inches="tight")
    fig.savefig(f"{out_file}.pdf", dpi=500, bbox_inches="tight", transparent=True)
    plt.close(fig)

    print(f"Saved TPM median plot to {out_file}.pdf")


def add_tpm_window_medians_to_df(
    df,
    window,
    tpm_log_col="log2_ref_tpm",
    out_col="median_window_refTPM",
):
    """
    Sliding-window medians per chromosome using only rows with finite TPM,
    but returns a column aligned to the FULL df (missing TPM rows remain NaN).
    """
    if "chrom" not in df.columns or "position" not in df.columns:
        raise KeyError("df must contain 'chrom' and 'position'.")
    if tpm_log_col not in df.columns:
        raise KeyError(f"df must contain '{tpm_log_col}'.")

    # full-length output
    out = np.full(len(df), np.nan, dtype=float)

    # only rows where TPM exists
    df_ok = df[np.isfinite(df[tpm_log_col])].copy()

    for chrom, sub in df_ok.groupby("chrom"):
        sub_sorted = sub.sort_values("position")
        idx = sub_sorted.index.to_numpy()

        vals = sub_sorted[tpm_log_col].to_numpy(dtype=float)

        med = sliding_median(vals, window)  # length n-window+1
        pad = len(vals) - len(med)
        med_full = np.concatenate([med, np.full(pad, np.nan)])

        out[idx] = med_full  # assign back using original indices

    df[out_col] = out
    return df

# --------------------------------------------------------------------
# Main
# --------------------------------------------------------------------
def main():
    args = get_args()

    print("Loading GFF...")
    gff = load_gtf(args.gff) # drop chrom to avoid issues during merge

    # Ensure one row per gene
    gff = (
        gff
        .sort_values(["gene", "chrom", "position"])
        .drop_duplicates(subset="gene", keep="first")
    )

    # because counts has chrom column too:
    gff = gff.drop(columns=["chrom"])

    # extract sample name from args.counts (/path/path{sample}_ase_by_reads_merged.txt)
    sample_name = args.counts.split("/")[-1].split("_ase_by_reads_merged.txt")[0]
    print(f"Sample name: {sample_name}")

    print("Loading counts...")
    counts = pd.read_csv(args.counts, sep="\t", comment="#")

    # Merge and sort by chromosome + position
    df = gff.merge(counts, on="gene", how="left")
    df = df.sort_values(["chrom", "position"])

    # -------- Using TPM --------

    # Compute log2 TPM
    df = compute_log2_tpm(df, tpm_col="ref_tpm", out_col="log2_ref_tpm", pseudocount=0.1)

    # (Optional) store window medians as a column, like your ASE pipeline
    df = add_tpm_window_medians_to_df(df, window=args.window, tpm_log_col="log2_ref_tpm",
                                    out_col="median_window_refTPM")

    # Plot
    plot_sliding_medians_by_chrom_tpm(
        df,
        window=args.window,
        tpm_log_col="log2_ref_tpm",
        out_file=f"{args.out}{sample_name}_medians_by_chromosomes_refTPM",
        y_lim=None  # or e.g. (0, 10) depending on your TPM scale
    )

    # ---------- Using ASE ---------

    # print("Computing ASE...")
    # df = compute_ase(df)

    # # Sliding stats per chromosome, using GLOBAL background
    # window = args.window
    # all_medians = []
    # all_pvals = []

    # # global ASE vector (used as background)
    # ase_all = df["ase"].values

    # for chrom, sub in df.groupby("chrom"):
    #     ase_vals = sub["ase"].values
    #     print(f"Processing chrom {chrom} with {len(sub)} genes...")
    #     med, pvals = sliding_stats(ase_vals, ase_all, window)

    #     # Pad with NaN so output has same number of rows as genes
    #     pad = len(sub) - len(med)
    #     med_full  = np.concatenate([med,   [np.nan] * pad])
    #     pval_full = np.concatenate([pvals, [np.nan] * pad])

    #     all_medians.append(med_full)
    #     all_pvals.append(pval_full)

    # df["median_window"] = np.concatenate(all_medians)
    # df["wilcoxon_p"]    = np.concatenate(all_pvals)
    # # add bonferroni correction column
    # n_ase = df["ase"].notna().sum()
    # if n_ase > 0:
    #     df["wilcoxon_p_bonferroni"] = df["wilcoxon_p"] * n_ase
    # else:
    #     df["wilcoxon_p_bonferroni"] = np.nan

    # explain_row(df, "chr1", "WASH7P", window=window)

    # print("Saving output...")
    # df.to_csv(f"{args.out}{sample_name}_aneuploidy_data.tsv", sep="\t", index=False)

    # print("Plotting medians by chromosome...")
    # plot_medians_by_chrom(df, window, out_file=f"{args.out}{sample_name}_medians_by_chromosomes")
    # print("Plotting median and p-values per chromosome...")
    # plot_median_and_p_per_chrom(df, window, out_prefix=f"{args.out}{sample_name}_chr_medians_pvals")

    # print("Done.")

if __name__ == "__main__":
    main()



# --------------------------------------------------------------------
# Load GFF2/GFF3 (simple parser assuming gene entries)
# --------------------------------------------------------------------
# def load_gff(gff_file):
#     rows = []
#     with open(gff_file) as f:
#         for line in f:
#             if line.startswith("#"):
#                 continue

#             parts = line.rstrip("\n").split("\t")
#             if len(parts) < 9:
#                 continue

#             chrom, source, feature, start, end, score, strand, phase, attributes = parts

#             # Only gene rows
#             if feature != "gene":
#                 continue

#             gene_name = None
#             gene_id = None

#             # Parse attributes
#             for field in attributes.split(";"):
#                 field = field.strip()
#                 if field.startswith("gene_name="):
#                     gene_name = field.split("=", 1)[1]
#                 elif field.startswith("gene_id="):
#                     gene_id = field.split("=", 1)[1]

#             # Fallback if gene_name missing
#             if gene_name is None:
#                 gene_name = gene_id

#             rows.append([gene_name, chrom, int(start)])

#     return pd.DataFrame(rows, columns=["gene", "chrom", "position"])
