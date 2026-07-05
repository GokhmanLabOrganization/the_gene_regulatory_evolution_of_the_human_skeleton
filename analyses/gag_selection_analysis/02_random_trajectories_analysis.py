import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from scipy.stats import uniform, kstest
import argparse
import sys
import os
from ast import literal_eval

# ==========================================
# 1. Unified Random Generation Logic
# ==========================================

def generate_random_fc_trajectories_matrix(df_pool, num_trajectories, n_per_trajectory, pool_col, min_nonzero=0, replace=False):
    """
    Generates trajectories directly into a NumPy Matrix (Float32).
    
    Args:
        replace (bool): 
            If False (Default) -> Sampling WITHOUT replacement (Original logic).
            If True -> Sampling WITH replacement (Optimized simple logic).
    """
    mode_str = "WITH REPLACEMENT" if replace else "WITHOUT REPLACEMENT"
    print(f"   > Optimization: Converting pool to Numpy for speed...")
    
    # Use Float32 for memory efficiency while maintaining sufficient precision
    all_values = df_pool[pool_col].to_numpy(dtype=np.float32) 
    total_pool_size = len(all_values)
    all_indices = np.arange(total_pool_size)
    
    nz_indices = np.flatnonzero(all_values != 0)
    
    # --- Validation ---
    if replace:
        # In replacement mode, we just need at least one non-zero if min_nonzero > 0
        if len(nz_indices) == 0 and min_nonzero > 0:
            raise ValueError("Error: Pool has no non-zero values, cannot satisfy min_nonzero.")
    else:
        # In NO-replacement mode, we need enough UNIQUE non-zeros
        if len(nz_indices) < min_nonzero:
            raise ValueError(
                f"Error: Pool has only {len(nz_indices)} non-zero values, "
                f"but you requested {min_nonzero} unique non-zeros per trajectory."
            )

    # --- Pre-allocation ---
    print(f"   > Generating {num_trajectories} trajectories ({mode_str})...")
    matrix = np.zeros((num_trajectories, n_per_trajectory), dtype=np.float32)

    for i in range(num_trajectories):
        
        # --- Step A: Head (Non-Zeros) ---
        # Pick indices for the required non-zero elements
        selected_nz_idx = np.random.choice(nz_indices, size=min_nonzero, replace=replace)
        
        # --- Step B: Tail (Remainder) ---
        needed = n_per_trajectory - min_nonzero
        
        if replace:
            # == WITH REPLACEMENT LOGIC ==
            # Simple: Sample from ALL indices. No need to mask anything.
            selected_remainder_idx = np.random.choice(all_indices, size=needed, replace=True)
            
            # Combine
            final_indices = np.concatenate([selected_nz_idx, selected_remainder_idx])
            # No shuffle needed for final_indices because order is already random from choice
            
        else:
            # == WITHOUT REPLACEMENT LOGIC ==
            # Complex: Must mask out the indices already picked in Step A
            # 1. Shuffle Head to randomize positions within the non-zero block (optional but good practice)
            np.random.shuffle(selected_nz_idx)
            
            # 2. Create Mask
            mask = np.ones(total_pool_size, dtype=bool)
            mask[selected_nz_idx] = False
            remaining_indices = all_indices[mask]
            
            # 3. Sample Tail from remaining
            if len(remaining_indices) < needed:
                 # Edge case: Take everything left
                 selected_remainder_idx = remaining_indices
            else:
                 selected_remainder_idx = np.random.choice(remaining_indices, size=needed, replace=False)
            
            np.random.shuffle(selected_remainder_idx)
            
            # 4. Combine
            final_indices = np.concatenate([selected_nz_idx, selected_remainder_idx])

        # --- Fill Matrix Row ---
        matrix[i, :] = all_values[final_indices]
        
    return matrix

# ==========================================
# 2. Statistical Analysis Logic (Fastest)
# ==========================================

def compute_last_position_pvalues_matrix(df_gene_pool, random_cumsum_matrix, id_col, traj_col, use_abs=True):
    """
    Calculates P-values using matrix slicing (Vectorized).
    Uses the "Background Matrix" approach for maximum speed.
    """
    # Pre-compute absolute values if needed (Vectorized on the whole matrix)
    # This trades RAM for Speed (faster lookups inside the loop)
    if use_abs:
        background_matrix = np.abs(random_cumsum_matrix)
    else:
        background_matrix = random_cumsum_matrix
        
    num_random_samples = background_matrix.shape[0]
    max_len = background_matrix.shape[1]

    results = []
    
    print("   > Calculating p-values (Vectorized Fast Mode)...")
    
    for idx, row in df_gene_pool.iterrows():
        gene_id = row.get(id_col, f"Gene_{idx}")
        fc = row.get(traj_col)
        
        if not isinstance(fc, (list, tuple, np.ndarray)) or len(fc) == 0:
            continue
        
        L = len(fc)
        
        if L > max_len: 
            continue
            
        last_val = fc[-1]
        if pd.isna(last_val): 
            continue

        target_val = abs(last_val) if use_abs else last_val
        
        # Access the specific column (Time step L-1) across all random trajectories
        # This is instantaneous with a matrix
        random_dist_at_L = background_matrix[:, L-1]
        
        # Vectorized comparison
        count_extreme = np.sum(random_dist_at_L >= target_val)
        
        p = (count_extreme + 1) / (num_random_samples + 1)
        
        results.append({
            id_col: gene_id,
            "last_value": last_val,
            "p_value": p,
            "L": L
        })

    return pd.DataFrame(results).sort_values("p_value").reset_index(drop=True)

# ==========================================
# 3. Plotting Functions 
# ==========================================

def plot_random_and_genes_matrix(random_matrix, genes_fc_lists, output_path):
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plot a subset of random background to avoid clutter
    max_plot = 500000
    subset = random_matrix[:max_plot, :]
    
    # Create X axis
    x_axis = np.arange(1, subset.shape[1] + 1)
    
    # Plot random background
    ax.plot(x_axis, subset.T, color='lightblue', linestyle='--', alpha=0.1, linewidth=0.5)

    # Plot genes
    for path in genes_fc_lists:
        if path and len(path) > 0:
            y = [x for x in path if isinstance(x, (int, float)) and not np.isnan(x)]
            ax.plot(np.arange(1, len(y)+1), y, color='purple', alpha=0.3, linewidth=0.5)

    ax.set_title(f"Random vs Gene Trajectories")
    ax.set_xlabel("Index (Sequence Position)")
    ax.set_ylabel("Cumulative logFC")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def plot_custom_qq(p_values, output_path):
    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    sm.qqplot(p_values, dist=uniform, line="45", ax=ax, marker='.', markersize=2, markeredgewidth=0, alpha=0.8)
    for line in ax.get_lines():
        line.set_linewidth(1.0)
    plt.title("QQ Plot of P-values")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def plot_log_qq(p_values, output_path):
    p_valid = p_values[~np.isnan(p_values)]
    p_valid = p_valid[(p_valid > 0) & (p_valid <= 1)]
    if len(p_valid) == 0: return

    ks_stat, ks_pval = kstest(p_valid, 'uniform')
    p_observed = np.sort(p_valid)
    n = len(p_observed)
    p_expected = np.arange(1, n + 1) / (n + 1)
    log_observed = -np.log10(p_observed)
    log_expected = -np.log10(p_expected)

    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    ax.plot(log_expected, log_observed, marker='.', linestyle='none', markersize=3, alpha=0.6, color='blue')
    max_val = max(log_expected.max(), log_observed.max())
    ax.plot([0, max_val], [0, max_val], color='red', linewidth=1.0, linestyle='--')
    
    stats_text = f"KS Test:\nD = {ks_stat:.4f}\np = {ks_pval:.2e}"
    ax.text(0.05, 0.85, stats_text, transform=ax.transAxes, fontsize=12, 
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.9))
            
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def plot_pvalue_distribution(csv_path, output_path=None, bins=50):
    if not os.path.exists(csv_path): return
    try: df = pd.read_csv(csv_path)
    except: return
    
    if 'p_value' not in df.columns: return
    p_values = df['p_value'].dropna()
    total_genes = len(p_values)
    if total_genes == 0: return

    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    counts, _, _ = ax.hist(p_values, bins=bins, color='skyblue', edgecolor='black', alpha=0.7)
    
    expected_count = total_genes / bins
    ax.axhline(expected_count, color='red', linestyle='--', linewidth=2, label='Expected (Uniform)')
    
    significant_05 = (p_values < 0.05).sum()
    stats_text = f"Total: {total_genes}\np<0.05: {significant_05}"
    
    ax.text(0.95, 0.95, stats_text, transform=ax.transAxes, fontsize=12, 
            verticalalignment='top', horizontalalignment='right', 
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.9))
            
    ax.set_title("P-value Distribution")
    ax.set_xlabel("P-value")
    ax.set_ylabel("Count")
    ax.legend()
    
    if output_path is None: 
        output_path = csv_path.replace(".csv", "_hist.png")
    
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

# ==========================================
# 4. Main Execution Block
# ==========================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run trajectory analysis (With or Without Replacement)")
    parser.add_argument('--genes_file', type=str, required=True)
    parser.add_argument('--pool_file', type=str, required=True)
    parser.add_argument('--output_dir', type=str, default="./results")
    parser.add_argument('--num_random', type=int, default=1000)
    parser.add_argument('--traj_len', type=int, default=None)
    parser.add_argument('--seed', type=int, default=None)
    parser.add_argument('--min_nonzero', type=int, default=0)
    parser.add_argument('--force_new', action='store_true')
    parser.add_argument('--id_col', type=str, default="NCBI_Gene_ID")
    parser.add_argument('--traj_col', type=str, default="logFC_active_scuffled_cum")
    parser.add_argument('--pool_col', type=str, default="logFC_derived_vs_ancestral")
    
    # NEW ARGUMENT: Flag for replacement
    # If present, args.with_replacement is True. If absent, it is False.
    parser.add_argument('--with_replacement', action='store_true', 
                        help="If set, random sampling is done WITH replacement. Default is WITHOUT replacement.")

    args = parser.parse_args()

    if args.seed is not None:
        np.random.seed(args.seed)

    os.makedirs(args.output_dir, exist_ok=True)
    
    # Determine mode for printing and filenames
    repl_mode = args.with_replacement
    mode_text = "WITH_REPLACEMENT" if repl_mode else "WITHOUT_REPLACEMENT"
    
    print(f"--- Starting Dynamic Analysis [{mode_text}] ---")
    
    # 1. Load Genes
    converters = {args.traj_col: literal_eval}
    try:
        df_gene_pool = pd.read_csv(args.genes_file, converters=converters)
    except Exception as e:
        sys.exit(f"Error loading genes file: {e}")

    if args.traj_len is None:
        print("Auto-detecting max trajectory length...")
        def get_len(x): return len(x) if isinstance(x, list) else 0
        args.traj_len = df_gene_pool[args.traj_col].apply(get_len).max()
        print(f"Detected max length: {args.traj_len}")

    # 2. Caching Logic
    # We include 'repl_True' or 'repl_False' in the filename to distinguish caches
    cache_filename = f"matrix_traj_N{args.num_random}_Len{args.traj_len}_Min{args.min_nonzero}_repl_{repl_mode}.npy"
    cache_path = os.path.join(args.output_dir, cache_filename)
    
    random_cumsum_matrix = None
    
    if os.path.exists(cache_path) and not args.force_new:
        print(f"Found cached matrix at: {cache_path}")
        try:
            random_cumsum_matrix = np.load(cache_path)
            print("Successfully loaded matrix.")
        except Exception as e:
            print(f"Error loading cache: {e}. Will regenerate.")
            
    if random_cumsum_matrix is None:
        print(f"Loading pool from: {args.pool_file}")
        try:
            df_pool = pd.read_csv(args.pool_file)
        except Exception as e:
            sys.exit(f"Error loading pool file: {e}")

        # A. Generate RAW Matrix
        raw_matrix = generate_random_fc_trajectories_matrix(
            df_pool, 
            args.num_random, 
            args.traj_len, 
            args.pool_col,
            min_nonzero=args.min_nonzero,
            replace=repl_mode  # Pass the flag here
        )
        
        # B. Vectorized Cumsum
        print("Calculating cumulative sums (Vectorized)...")
        np.cumsum(raw_matrix, axis=1, out=raw_matrix)
        random_cumsum_matrix = raw_matrix
        
        # Save
        print(f"Saving generated matrix to: {cache_path}")
        np.save(cache_path, random_cumsum_matrix)
        
        # Cleanup
        del raw_matrix
        import gc
        gc.collect()

    # 3. Compute Statistics
    df_extreme = compute_last_position_pvalues_matrix(
        df_gene_pool,
        random_cumsum_matrix,
        id_col=args.id_col,
        traj_col=args.traj_col,
        use_abs=True
    )

    if df_extreme.empty:
        sys.exit("No valid genes found.")

    # 4. FDR Correction
    pvals = df_extreme["p_value"].values
    reject, pvals_fdr, _, _ = multipletests(pvals, alpha=0.05, method="fdr_bh")
    df_extreme["p_adj"] = pvals_fdr
    df_extreme["reject_fdr_0.05"] = reject

    # 5. Save Results
    output_csv = os.path.join(args.output_dir, "analysis_results.csv")
    df_extreme.to_csv(output_csv, index=False)
    print(f"Results saved to: {output_csv}")

    # 6. Plotting
    print("Generating plots...")
    plot_path_traj = os.path.join(args.output_dir, "trajectories_plot.png")
    plot_random_and_genes_matrix(
        random_cumsum_matrix, 
        df_gene_pool[args.traj_col].tolist(), 
        plot_path_traj
    )
    
    plot_path_qq_lin = os.path.join(args.output_dir, "qq_plot_linear.png")
    plot_custom_qq(pvals, plot_path_qq_lin)

    plot_path_qq_log = os.path.join(args.output_dir, "qq_plot_log.png")
    plot_log_qq(pvals, plot_path_qq_log)

    plot_path_hist = os.path.join(args.output_dir, "pvalue_histogram.png")
    plot_pvalue_distribution(output_csv, plot_path_hist)

    print("Done.")