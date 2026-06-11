#!/usr/bin/env python3
"""
Compute ref_TPM, alt_TPM and total_TPM from a counts file and a GFF file.

Inputs
------
1. Counts file (TSV) with at least columns:
   gene, chrom, ref_counts, alt_counts, no_ase_counts, ambig_ase_counts
   - May contain extra columns (ignored for TPM).
   - Lines starting with '#' are treated as comments and skipped.

2. GFF file path:
   - Must contain features of type "gene".
   - Gene length is end - start + 1.
   - Gene ID is taken from GFF attributes: 'ID=' or 'gene_id='.

3. Excluded regions:
   - A comma-separated string, e.g.:
       "chrM,chrUn,chr1:100000-200000,chr2:50000-60000"
   - Items without ':' are whole chromosomes to remove.
   - Items with ':' are intervals [start, end] (1-based, inclusive)
     on a chromosome.

4. Outpath:
   - Path to write the filtered counts + TPM columns (TSV).

Output
------
A TSV file with all original counts columns plus:
    ref_TPM, alt_TPM, total_TPM
for genes NOT overlapping excluded regions and with known length.
"""

import argparse
import sys
import pandas as pd
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute ref_TPM, alt_TPM and total_TPM from counts and GFF."
    )
    parser.add_argument("--counts", required=True,
                        help="Counts file path (TSV; lines starting with # ignored).")
    parser.add_argument("--gff", required=True,
                        help="GFF file with gene annotations.")
    parser.add_argument("--exclude", default="",
                        help=(
                            "Chromosomes/regions to remove, comma-separated. "
                            "Examples: 'chrM,chrUn,chr1:100000-200000'. "
                            "Items without ':' are whole chromosomes."
                        ))
    parser.add_argument("--out", required=True,
                        help="Output TSV path for filtered counts with TPM columns.")
    return parser.parse_args()


def parse_excluded_regions(spec):
    """
    Parse exclusion specification string.

    Returns:
        full_chroms: set of chromosome names to exclude entirely.
        regions: dict of chrom -> list of (start, end) intervals to exclude.
    """
    full_chroms = set()
    regions = {}

    spec = spec.strip()
    if not spec:
        return full_chroms, regions

    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        if ":" in token:
            chrom, rest = token.split(":", 1)
            try:
                start_str, end_str = rest.split("-", 1)
                start = int(start_str.replace(",", ""))
                end = int(end_str.replace(",", ""))
            except ValueError:
                raise ValueError(f"Cannot parse region '{token}'. Expected format 'chr:start-end'.")
            if chrom not in regions:
                regions[chrom] = []
            regions[chrom].append((start, end))
        else:
            full_chroms.add(token)

    print(f"Excluding chromosomes: {full_chroms}", file=sys.stderr)
    print(f"Excluding regions: {regions}", file=sys.stderr)
    return full_chroms, regions


def parse_gff_attributes(attr_str):
    """
    Parse GFF3 attribute column into a dict.
    Handles simple 'key=value;key2=value2' and also 'key "value"' cases.

    Returns:
        dict of {attribute_key: attribute_value}
    """
    attrs = {}
    if not attr_str or attr_str == ".":
        return attrs

    for item in attr_str.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, val = item.split("=", 1)
            attrs[key.strip()] = val.strip()
        else:
            # GTF-style: key "value"
            parts = item.split()
            if len(parts) >= 2:
                key = parts[0]
                val = parts[1].strip().strip('"')
                attrs[key] = val
    return attrs


# def load_gene_lengths_from_gtf_by_symbol(gtf_path, name_attr="gene_name"):
#     """
#     Read GTF/GFF and extract gene lengths and coordinates, keyed by gene symbol.

#     Returns:
#         gene_info: DataFrame with columns
#           gene_name, gff_chrom, gff_start, gff_end, gene_length
#     """
#     gene_records = []

#     with open(gtf_path, "r") as fh:
#         for line in fh:
#             if not line.strip() or line.startswith("#"):
#                 continue

#             parts = line.rstrip("\n").split("\t")
#             if len(parts) < 9:
#                 continue

#             chrom, source, feature_type, start, end, score, strand, phase, attrs_str = parts

#             # Only "gene" features
#             if feature_type != "gene":
#                 continue

#             attrs = parse_gff_attributes(attrs_str)
#             gene_name = attrs.get(name_attr)
#             if gene_name is None:
#                 continue

#             try:
#                 start_i = int(start)
#                 end_i = int(end)
#             except ValueError:
#                 continue

#             length = end_i - start_i + 1
#             if length <= 0:
#                 continue

#             gene_records.append(
#                 {
#                     "gene_name": gene_name,
#                     "gff_chrom": chrom,
#                     "gff_start": start_i,
#                     "gff_end": end_i,
#                     "gene_length": length,
#                 }
#             )

#     if not gene_records:
#         raise RuntimeError("No gene features with gene_name found in GTF/GFF.")

#     gene_info = pd.DataFrame(gene_records)

#     # If a symbol appears multiple times (rare), keep the longest one
#     gene_info = (
#         gene_info.sort_values(["gene_name", "gene_length"], ascending=[True, False])
#         .drop_duplicates(subset="gene_name", keep="first")
#         .reset_index(drop=True)
#     )

#     return gene_info

def load_gene_lengths_from_gtf_exons_by_symbol(gtf_path, name_attr="gene_name"):
    """
    Build gene lengths from exon features in a GTF/GFF, keyed by gene symbol.

    For each gene_name:
      - collect all exon intervals per chromosome
      - merge overlapping/adjacent intervals
      - gene_length = sum of merged exon lengths (union)
      - gff_start = min exon start
      - gff_end   = max exon end

    Returns:
        gene_info: DataFrame with columns:
          gene_name, gff_chrom, gff_start, gff_end, gene_length
    """
    exons_by_gene = defaultdict(lambda: defaultdict(list))

    with open(gtf_path, "r") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            chrom, source, feature_type, start, end, score, strand, phase, attrs_str = parts

            if feature_type != "exon":
                continue

            attrs = parse_gff_attributes(attrs_str)
            gene_name = attrs.get(name_attr)
            if gene_name is None:
                continue

            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue
            if end_i < start_i:
                continue

            exons_by_gene[gene_name][chrom].append((start_i, end_i))

    gene_records = []

    for gene_name, chrom_dict in exons_by_gene.items():
        for chrom, intervals in chrom_dict.items():
            if not intervals:
                continue

            intervals_sorted = sorted(intervals, key=lambda x: x[0])

            merged_intervals = []
            cur_start, cur_end = intervals_sorted[0]

            for s, e in intervals_sorted[1:]:
                if s <= cur_end + 1:
                    # overlap or directly adjacent
                    cur_end = max(cur_end, e)
                else:
                    merged_intervals.append((cur_start, cur_end))
                    cur_start, cur_end = s, e
            merged_intervals.append((cur_start, cur_end))

            union_len = sum(e - s + 1 for s, e in merged_intervals)
            gene_start = min(s for s, e in merged_intervals)
            gene_end   = max(e for s, e in merged_intervals)

            gene_records.append(
                {
                    "gene_name": gene_name,
                    "gff_chrom": chrom,
                    "gff_start": gene_start,
                    "gff_end": gene_end,
                    "gene_length": union_len,
                }
            )

    if not gene_records:
        raise RuntimeError("No exon-based gene records could be built from GTF.")

    gene_info = pd.DataFrame(gene_records)

    # if same symbol appears multiple times, keep longest entry
    gene_info = (
        gene_info.sort_values(["gene_name", "gene_length"], ascending=[True, False])
        .drop_duplicates(subset="gene_name", keep="first")
        .reset_index(drop=True)
    )

    return gene_info



def load_counts(counts_path):
    """
    Load only the needed columns from counts file, and also return header
    comment lines (starting with '#') so we can keep them in the output.

    Returns:
        df, header_comments
    """
    needed_cols = [
        "gene",
        "chrom",
        "ref_counts",
        "alt_counts",
        "no_ase_counts",
        "ambig_ase_counts",
    ]

    # collect comment/header lines
    header_comments = []
    with open(counts_path, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                header_comments.append(line.rstrip("\n"))
            else:
                # first non-# line -> done
                break

    # now read the table; `comment='#'` will ignore those lines anyway
    df = pd.read_csv(
        counts_path,
        sep="\t",
        comment="#",
        usecols=lambda c: c in needed_cols,
        dtype=str,
    )

    # numeric columns
    for c in ["ref_counts", "alt_counts", "no_ase_counts", "ambig_ase_counts"]:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(float)

    return df, header_comments


def apply_exclusions(df, full_chroms, regions):
    """
    Remove genes based on excluded chromosomes and regions.

    df must have:
      gff_chrom, gff_start, gff_end

    Returns filtered DataFrame (copy).
    """
    mask = pd.Series(True, index=df.index)

    # Exclude entire chromosomes
    if full_chroms:
        mask &= ~df["gff_chrom"].isin(full_chroms)

    # Exclude regions by overlap
    for chrom, intervals in regions.items():
        # Only rows on this chromosome
        on_chrom = df["gff_chrom"].eq(chrom)
        if not on_chrom.any():
            continue
        # Build a per-row mask for this chromosome
        start = df.loc[on_chrom, "gff_start"]
        end = df.loc[on_chrom, "gff_end"]
        # Initially, nothing excluded on this chrom
        exclude_here = pd.Series(False, index=start.index)

        for (s, e) in intervals:
            # Overlap condition: gene_end >= s and gene_start <= e
            overlap = (end >= s) & (start <= e)
            exclude_here |= overlap

        # Combine into global mask
        mask.loc[exclude_here.index] &= ~exclude_here

    filtered = df.loc[mask].copy()
    return filtered


def compute_tpm(df):
    """
    Compute ref_TPM, alt_TPM, total_TPM given:

    df columns:
        ref_counts, alt_counts, non_ase_counts, ambig_counts, gene_length

    Strategy:
        length_kb = gene_length / 1000
        total_counts = ref + alt
        RPK_total = total_counts / length_kb
        scaling_factor = sum(RPK_total) / 1e6
        total_TPM = RPK_total / sum(RPK_total) * 1e6
        ref_TPM = ref_counts / length_kb / sum(RPK_total) * 1e6
        alt_TPM = alt_counts / length_kb / sum(RPK_total) * 1e6

    So ref_TPM + alt_TPM == total_TPM (up to float rounding).
    """
    length_kb = df["gene_length"].astype(float) / 1000.0
    # Avoid division by zero
    length_kb[length_kb <= 0] = float("nan")

    ref_counts = df["ref_counts"].astype(float)
    alt_counts = df["alt_counts"].astype(float)
    non_ase_counts = df["no_ase_counts"].astype(float)
    ambig_ase_counts = df["ambig_ase_counts"].astype(float)
    total_counts = ref_counts + alt_counts + non_ase_counts + ambig_ase_counts

    RPK_total = total_counts / length_kb

    denom = RPK_total.sum()
    if denom <= 0 or pd.isna(denom):
        raise RuntimeError("Total RPK is zero or NaN; cannot compute TPM.")

    factor = 1e6 / denom

    df["total_TPM"] = RPK_total * factor
    df["ref_TPM"] = (ref_counts / length_kb) * factor
    df["alt_TPM"] = (alt_counts / length_kb) * factor

    return df


def main():
    args = parse_args()

    full_chroms, regions = parse_excluded_regions(args.exclude)

    # print(f"Loading gene lengths from GTF: {args.gff}", file=sys.stderr)
    # gene_info = load_gene_lengths_from_gtf_by_symbol(args.gff)

    print(f"Loading gene lengths from GTF (exon-based): {args.gff}", file=sys.stderr)
    gene_info = load_gene_lengths_from_gtf_exons_by_symbol(args.gff)

    print(f"Loading counts from: {args.counts}", file=sys.stderr)
    counts, header_comments = load_counts(args.counts) 

    # Merge on gene symbol
    merged = counts.merge(
        gene_info,
        left_on="gene",
        right_on="gene_name",
        how="left",
        validate="m:1",
    )

    # print the missing gene names 
    missing_mask = merged["gene_length"].isna()
    if missing_mask.any():
        missing_genes = merged.loc[missing_mask, "gene"].unique()
        print(
            f"Warning: {len(missing_genes)} genes from counts have no gene_length in GTF and will be dropped:",
            file=sys.stderr,
        )
        # print them as a comma-separated list (or one per line if you prefer)
        print("  " + ", ".join(sorted(map(str, missing_genes))), file=sys.stderr)

        # now drop them
        merged = merged[~missing_mask].copy()

    # # Alternatively, just print count of missing genes and drop them
    # missing_len = merged["gene_length"].isna().sum()
    # if missing_len > 0:
    #     print(
    #         f"Warning: {missing_len} genes from counts have no gene_length in GTF and will be dropped.",
    #         file=sys.stderr,
    #     )
    #     merged = merged[merged["gene_length"].notna()].copy()

    merged = merged.drop(columns=["gene_name"])

    if full_chroms or regions:
        print("Applying chromosome/region exclusions...", file=sys.stderr)
        before = len(merged)
        merged = apply_exclusions(merged, full_chroms, regions)
        after = len(merged)
        print(f"Excluded {before - after} genes based on provided regions.", file=sys.stderr)

    print("Computing TPM...", file=sys.stderr)
    merged = compute_tpm(merged)


    # drop GFF helper columns if not needed in output
    merged = merged.drop(columns=["gff_chrom", "gff_start", "gff_end", "gene_length"])

    print(f"Writing output to: {args.out}", file=sys.stderr)
    with open(args.out, "w") as out_f:
        # write the original header comment lines
        for line in header_comments:
            out_f.write(line + "\n")
        # then write the table into the SAME handle
        merged.to_csv(out_f, sep="\t", header=True, index=False)


if __name__ == "__main__":
    main()
