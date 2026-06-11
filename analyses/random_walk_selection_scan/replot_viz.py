import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
from ast import literal_eval
import os
import sys
import argparse
from matplotlib.lines import Line2D

# ==========================================
# 1. Visualization Logic
# ==========================================

def plot_random_and_genes_matrix(random_matrix, df_genes, traj_col, sig_col, output_path):
    """
    Plots the random background, statistical thresholds, and gene trajectories.
    
    Updates:
    - Optimized for ~3500 genes (increased alpha and linewidth).
    - Uses 'dodgerblue' for non-significant genes for better contrast.
    """
    print(f"   > Initializing plot canvas...")
    
    # Create figure with high DPI for publication quality
    fig, ax = plt.subplots(figsize=(12, 8), dpi=300)
    
    # --- A. Random Background (Subset) ---
    # We plot only a subset of random trajectories to avoid memory crashes 
    # and solid blocks of color.
    max_plot_random = 100000
    limit = min(max_plot_random, random_matrix.shape[0])
    subset = random_matrix[:limit, :]
    x_axis_rand = np.arange(1, subset.shape[1] + 1)
    
    print(f"   > Plotting background ({limit} random trajectories)...")
    # Light grey background
    ax.plot(x_axis_rand, subset.T, color='lightgrey', linestyle='-', alpha=0.05, linewidth=0.5, zorder=1)

    # --- B. Statistical Thresholds (Full Data) ---
    # Critical: We calculate the threshold using the FULL matrix (e.g. 500k) 
    # to ensure statistical accuracy, even if we only plotted a subset.
    print("   > Calculating 95th percentile thresholds (P=0.05)...")
    abs_matrix = np.abs(random_matrix)
    
    # Calculate quantile per time-step (column)
    threshold_95 = np.percentile(abs_matrix, 95, axis=0) 
    
    # Plot Positive and Negative bounds (Red Dashed Line)
    ax.plot(x_axis_rand, threshold_95, color='red', linestyle='--', linewidth=1.2, alpha=0.9, label='P=0.05 Threshold', zorder=5)
    ax.plot(x_axis_rand, -threshold_95, color='red', linestyle='--', linewidth=1.2, alpha=0.9, zorder=5)

    # --- C. Gene Trajectories ---
    print(f"   > Plotting {len(df_genes)} gene trajectories...")
    
    # --- STYLING CONFIGURATION ---
    # Updated for better visibility with ~3500 genes:
    # 1. Color: 'dodgerblue' (brighter than cornflowerblue).
    # 2. Alpha: 0.2 (more opaque).
    # 3. Linewidth: 1.0 (thicker).
    style_nonsig = {'color': 'dodgerblue', 'alpha': 0.2, 'linewidth': 1.0, 'zorder': 2}
    
    # Significant genes: Darker, thicker, and more on top.
    style_sig =    {'color': 'navy',       'alpha': 0.8, 'linewidth': 1.2, 'zorder': 3}
    
    count_sig = 0
    count_nonsig = 0

    for _, row in df_genes.iterrows():
        path = row[traj_col]
        is_sig = row[sig_col]
        
        # Validation
        if not isinstance(path, (list, np.ndarray)) or len(path) == 0:
            continue
            
        y = np.array(path)
        x = np.arange(1, len(y) + 1)
        
        if is_sig:
            ax.plot(x, y, **style_sig)
            count_sig += 1
        else:
            ax.plot(x, y, **style_nonsig)
            count_nonsig += 1

    # --- D. Final Touches & Legend ---
    print("   > Finalizing layout...")
    
    # Custom Legend Handles
    legend_elements = [
        Line2D([0], [0], color='lightgrey', lw=2, label='Random Model'),
        Line2D([0], [0], color='dodgerblue', lw=2, label=f'Non-Sig (FDR>0.05): {count_nonsig}'),
        Line2D([0], [0], color='navy', lw=2, label=f'Significant (FDR<0.05): {count_sig}'),
        Line2D([0], [0], color='red', linestyle='--', lw=2, label='P=0.05 Threshold')
    ]

    ax.set_title("Evolutionary Trajectories: Genes vs Random Model", fontsize=14)
    ax.set_xlabel("Sequence Position", fontsize=12)
    ax.set_ylabel("Cumulative logFC", fontsize=12)
    ax.legend(handles=legend_elements, loc='upper left', fontsize=10)
    ax.grid(True, which='both', linestyle=':', linewidth=0.5, alpha=0.3)
    
    # Auto-adjust layout to prevent cutting off labels
    plt.tight_layout()
    
    # Save
    plt.savefig(output_path)
    print(f"   > Done. Plot saved to: {output_path}")
    plt.close()

# ==========================================
# 2. Argument Parsing & Execution
# ==========================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Re-plot visualization from cached data.")
    
    # Required Arguments
    parser.add_argument('--genes_file', type=str, required=True, help="Path to original genes CSV file")
    parser.add_argument('--results_file', type=str, required=True, help="Path to analysis_results.csv (with p-values)")
    parser.add_argument('--matrix_file', type=str, required=True, help="Path to the .npy random matrix file")
    
    # Optional Arguments (with defaults)
    parser.add_argument('--output_plot', type=str, default="replot_viz.png", help="Output path for the image")
    parser.add_argument('--id_col', type=str, default="NCBI_Gene_ID", help="Column name for Gene ID")
    parser.add_argument('--traj_col', type=str, default="logFC_active_scuffled_cum", help="Column name for trajectory list")
    parser.add_argument('--sig_col', type=str, default="reject_fdr_0.05", help="Column name for significance boolean")

    args = parser.parse_args()
    
    # --- Validations ---
    if not os.path.exists(args.genes_file):
        sys.exit(f"Error: Genes file not found at: {args.genes_file}")
    if not os.path.exists(args.results_file):
        sys.exit(f"Error: Results file not found at: {args.results_file}")
    if not os.path.exists(args.matrix_file):
        sys.exit(f"Error: Matrix file not found at: {args.matrix_file}")

    # --- Loading Data ---
    print(f"1. Loading random matrix: {args.matrix_file}")
    try:
        random_matrix = np.load(args.matrix_file)
    except Exception as e:
        sys.exit(f"Error loading numpy file: {e}")

    print(f"2. Loading genes: {args.genes_file}")
    # Using literal_eval to parse the trajectory string "[0.1, 0.2]" into a python list
    df_genes_raw = pd.read_csv(args.genes_file, converters={args.traj_col: literal_eval})

    print(f"3. Loading results: {args.results_file}")
    df_results = pd.read_csv(args.results_file)

    # --- Merging ---
    print("4. Merging data...")
    # Ensure ID columns are treated as strings to avoid type mismatch (int vs str)
    if args.id_col in df_genes_raw.columns:
        df_genes_raw[args.id_col] = df_genes_raw[args.id_col].astype(str)
    
    if args.id_col in df_results.columns:
        df_results[args.id_col] = df_results[args.id_col].astype(str)

    # Merge: Keep only genes that exist in both files (Inner Join)
    df_plot = pd.merge(
        df_genes_raw[[args.id_col, args.traj_col]], 
        df_results[[args.id_col, args.sig_col]], 
        on=args.id_col, 
        how="inner"
    )

    if df_plot.empty:
        sys.exit("Error: Merged DataFrame is empty. Please check if Gene IDs match between the input file and the results file.")
    
    print(f"   > Data ready. Proceeding to plot {len(df_plot)} genes.")

    # --- Plotting ---
    plot_random_and_genes_matrix(
        random_matrix=random_matrix,
        df_genes=df_plot,
        traj_col=args.traj_col,
        sig_col=args.sig_col,
        output_path=args.output_plot
    )