# plot_utils.py

import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy.stats import uniform
import matplotlib.patches as mpatches
def describe_series(vals):
    """
    Return basic distribution statistics for a 1D pandas Series:
    count, min, max, mean, median, std.
    """
    return {
        "total": len(vals),
        "min": vals.min(),
        "max": vals.max(),
        "mean": vals.mean(),
        "median": vals.median(),
        "std": vals.std(),
    }


def build_summary_text(
    stats_all,
    stats_diff,
    header_all="Active",
    header_diff="Diff active",
    count_label_all="Total",
    count_label_diff="Total",
    value_label="",
):
    """
    Build a two-block text summary for the legend box on the central plot.

    Parameters
    ----------
    stats_all, stats_diff : dict
        Output of `describe_series` for the two distributions.
    header_all, header_diff : str
        Titles for the two blocks ("Active", "Diff active", etc.).
    count_label_all, count_label_diff : str
        How to label the count ("Total", "Genes", etc.).
    value_label : str
        Optional label appended after 'Min', 'Max', 'Mean', etc.
        e.g. value_label="mean" -> 'Min mean', 'Mean mean', etc.
        If empty, we just use 'Min', 'Max', 'Mean', ...
    """
    suffix = f" {value_label}" if value_label else ""

    text_all = (
        f"{header_all}:\n"
        f"  {count_label_all}: {stats_all['total']}\n"
        f"  Min{suffix}: {stats_all['min']:.2f}\n"
        f"  Max{suffix}: {stats_all['max']:.2f}\n"
        f"  Mean{suffix}: {stats_all['mean']:.3f}\n"
        f"  Median{suffix}: {stats_all['median']:.3f}\n"
        f"  Std Dev{suffix}: {stats_all['std']:.3f}\n"
    )

    text_diff = (
        f"{header_diff}:\n"
        f"  {count_label_diff}: {stats_diff['total']}\n"
        f"  Min{suffix}: {stats_diff['min']:.2f}\n"
        f"  Max{suffix}: {stats_diff['max']:.2f}\n"
        f"  Mean{suffix}: {stats_diff['mean']:.3f}\n"
        f"  Median{suffix}: {stats_diff['median']:.3f}\n"
        f"  Std Dev{suffix}: {stats_diff['std']:.3f}"
    )

    return text_all + "\n\n" + text_diff


def plot_center_tails(
    vals_all,
    vals_diff,
    xlabel,
    ylabel_center,
    ylabel_tails,
    label_all_center="Active",
    label_diff_center="Diff active",
    label_all_tails=None,
    label_diff_tails=None,
    center_title="Center (1–99% quantiles)",
    tails_title="Tails",
    summary_kwargs=None,
    print_prefix="Central range",
    lower_q=1,
    upper_q=99,
    bins_center_count=40,
    bins_tails_count=40,
    figsize=(16, 6),
    width_ratios=(2, 1),
):
    """
    Plot a 1x2 figure with:
    - left panel: central quantile range (overlapping histograms of vals_all / vals_diff)
    - right panel: tails (everything outside the central range)

    Parameters
    ----------
    vals_all, vals_diff : pandas Series
        Full and "diff" distributions to compare (NaNs are dropped inside).
    xlabel : str
        X-axis label for both panels.
    ylabel_center, ylabel_tails : str
        Y-axis labels for center and tails panels, respectively.
    label_all_center, label_diff_center : str
        Legend labels for the central panel.
    label_all_tails, label_diff_tails : str or None
        Legend labels for the tails panel. If None, defaults to
        "<label_all_center> (tails)" and "<label_diff_center> (tails)".
    center_title, tails_title : str
        Titles for the two panels.
    summary_kwargs : dict or None
        Keyword arguments forwarded to `build_summary_text` to control
        the summary box (headers, count labels, value_label, etc.).
    print_prefix : str
        Prefix for the printed description of the central range in stdout.
    lower_q, upper_q : float
        Quantiles used to define the "central" range.
    bins_center_count, bins_tails_count : int
        Number of bins for the central and tails histograms.
    figsize : tuple
        Figure size passed to plt.subplots.
    width_ratios : tuple
        Relative widths of the center and tails panels.
    """

    # Drop NaNs (we want the same behavior as your original code)
    vals_all = vals_all.dropna()
    vals_diff = vals_diff.dropna()

    # Summary statistics
    stats_all = describe_series(vals_all)
    stats_diff = describe_series(vals_diff)

    # Central quantile range
    lower, upper = np.percentile(vals_all, [lower_q, upper_q])

    central_all = vals_all[(vals_all >= lower) & (vals_all <= upper)]
    central_diff = vals_diff[(vals_diff >= lower) & (vals_diff <= upper)]

    tails_all = vals_all[(vals_all < lower) | (vals_all > upper)]
    tails_diff = vals_diff[(vals_diff < lower) | (vals_diff > upper)]

    # Console output (keeps your original debug info, just parameterized)
    print(f"{print_prefix}: [{lower:.2f}, {upper:.2f}]")
    print(f"Central all: {len(central_all)}, tails all: {len(tails_all)}")
    print(f"Central diff: {len(central_diff)}, tails diff: {len(tails_diff)}")

    # Bins
    bins_center = np.linspace(lower, upper, bins_center_count)
    tails_min = vals_all.min()
    tails_max = vals_all.max()
    bins_tails = np.linspace(tails_min, tails_max, bins_tails_count)

    # Figure + axes
    fig, (ax_center, ax_tails) = plt.subplots(
        1, 2, figsize=figsize, gridspec_kw={"width_ratios": width_ratios}
    )

    # ----- LEFT: central zoom -----
    ax_center.hist(
        central_all,
        bins=bins_center,
        edgecolor="black",
        alpha=0.5,
        color="skyblue",
        label=label_all_center,
    )

    ax_center.hist(
        central_diff,
        bins=bins_center,
        edgecolor="black",
        alpha=0.5,
        color="red",
        label=label_diff_center,
    )

    # Vertical mean lines
    ax_center.axvline(
        stats_all["mean"],
        color="blue",
        linestyle="--",
        linewidth=1.5,
        label=f"Mean {label_all_center} = {stats_all['mean']:.2f}",
    )
    ax_center.axvline(
        stats_diff["mean"],
        color="red",
        linestyle="--",
        linewidth=1.5,
        label=f"Mean {label_diff_center} = {stats_diff['mean']:.2f}",
    )

    ax_center.set_xticks(
        np.arange(np.floor(lower), np.ceil(upper) + 0.5, 0.5)
    )
    ax_center.set_xlabel(xlabel)
    ax_center.set_ylabel(ylabel_center)
    ax_center.set_title(center_title)
    ax_center.legend()

    # Summary text box (right side of central panel)
    if summary_kwargs is None:
        summary_kwargs = {}

    summary_text = build_summary_text(
        stats_all,
        stats_diff,
        **summary_kwargs,
    )

    ax_center.text(
        1.02,
        0.95,
        summary_text,
        transform= ax_center.transAxes,
        fontsize=11,
        va="top",
        fontweight="bold",
    )

    # ----- RIGHT: tails -----
    if label_all_tails is None:
        label_all_tails = f"{label_all_center} (tails)"
    if label_diff_tails is None:
        label_diff_tails = f"{label_diff_center} (tails)"

    ax_tails.hist(
        tails_all,
        bins=bins_tails,
        edgecolor="black",
        alpha=0.5,
        color="skyblue",
        label=label_all_tails,
    )

    ax_tails.hist(
        tails_diff,
        bins=bins_tails,
        edgecolor="black",
        alpha=0.5,
        color="red",
        label=label_diff_tails,
    )

    ax_tails.set_xlabel(xlabel)
    ax_tails.set_ylabel(ylabel_tails)
    ax_tails.set_title(tails_title)

    # Mark the central cutoffs
    ax_tails.axvline(lower, color="grey", linestyle="--", linewidth=1)
    ax_tails.axvline(upper, color="grey", linestyle="--", linewidth=1)

    ax_tails.legend()

    plt.tight_layout()
    plt.show()

    # In case you ever want these thresholds/stats programmatically
    return {
        "lower": lower,
        "upper": upper,
        "stats_all": stats_all,
        "stats_diff": stats_diff,
    }


def plot_enhancer_stack_by_pct(
vals,
    pct_bins,
    xlabel="Number of diff-active enhancers per gene",
    title="Distribution of diff-active enhancers per gene",
):
    """
    Plot a stacked bar chart with grouped colors for continuous visual blocks
    and a custom simplified legend.
    """

    # 1. Define Colors
    color_null = '#E8F5E9'   # Very Light Green
    color_mid  = '#66BB6A'   # Medium Green
    color_ext  = '#1B5E20'   # Dark Green

    # 2. Define Plotting Groups
    # We group specific bins together so they are plotted as a single block.
    # Note: We separate 'Extreme (Bottom)' and 'Extreme (Top)' so the stacking order 
    # remains correct (Null -> Low Extremes -> Middle -> High Extremes).
    plotting_groups = [
        {'label': 'Null',             'bins': ['null'],               'color': color_null},
        {'label': 'Extreme (Bottom)', 'bins': ['0%', '1–25%'],        'color': color_ext},
        {'label': 'Middle',           'bins': ['25–50%', '50–75%'],   'color': color_mid},
        {'label': 'Extreme (Top)',    'bins': ['75–99%', '100%'],     'color': color_ext},
    ]

    # Flatten list to ensure we calculate raw data for all expected bins
    all_possible_bins = [b for group in plotting_groups for b in group['bins']]

    # --- Data Preparation ---
    mask_valid = vals.notna() & pct_bins.notna()
    vals = vals[mask_valid]
    pct_bins = pct_bins[mask_valid]

    enhancer_counts = sorted(
        set(list(range(0, int(vals.dropna().max()) + 1)) + list(vals.dropna().unique()))
    )

    # Calculate raw counts for every fine-grained bin first
    data_raw = {pct: [] for pct in all_possible_bins}
    for enh_count in enhancer_counts:
        mask_count = vals == enh_count
        for pct in all_possible_bins:
            mask_pct = pct_bins == pct
            count = (mask_count & mask_pct).sum()
            data_raw[pct].append(count)

    # --- Plotting ---
    fig, ax = plt.subplots(figsize=(18, 6))
    x_pos = np.arange(len(enhancer_counts))
    bottom = np.zeros(len(enhancer_counts), dtype=int)

    # Loop through the logical GROUPS, not the individual bins
    for group in plotting_groups:
        # Sum heights of all bins belonging to this group
        group_heights = np.zeros(len(enhancer_counts), dtype=int)
        for bin_name in group['bins']:
            group_heights += np.array(data_raw[bin_name])
        
        # Plot the consolidated block
        ax.bar(
            x_pos,
            group_heights,
            color=group['color'],
            edgecolor='black', # Border draws around the whole group, not internal bins
            linewidth=0.8,
            bottom=bottom,
        )
        bottom += group_heights

    # --- Annotations & Styling ---
    # Total counts on top
    for i, total in enumerate(bottom):
        if total == 0: continue
        ax.text(x_pos[i], total + 0.5, str(int(total)), ha='center', va='bottom', fontsize=9)

    # Statistics box
    count_non_na = int(vals.dropna().shape[0])
    mean_enh = vals.mean()
    median_enh = vals.median()
    std_enh = vals.std()

    summary_text = (
        f"Genes with valid enhancer count: {count_non_na}\n\n"
        f"Mean enhancers per gene: {mean_enh:.2f}\n\n"
        f"Median enhancers per gene: {median_enh}\n\n"
        f"Std Dev enhancers per gene: {std_enh:.2f}"
    )
    ax.text(1.02, 0.95, summary_text, transform=ax.transAxes, fontsize=12, va='top', fontweight='bold')

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Number of genes")
    ax.set_title(title)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(enhancer_counts)

    # --- Custom Legend ---
    # Create manual handles to represent the 3 distinct color meanings
    legend_patches = [
        mpatches.Patch(facecolor=color_null, edgecolor='black', label='Null (NaN)'),
        mpatches.Patch(facecolor=color_mid,  edgecolor='black', label='Middle (25%–75%)'),
        mpatches.Patch(facecolor=color_ext,  edgecolor='black', label='Extreme (<25% or >75%)'),
    ]
    
    ax.legend(handles=legend_patches, title="Status / % Positive", loc='upper right')

    plt.tight_layout()
    return fig, ax

def plot_custom_qq(p_values, dist=uniform, title="QQ plot of p-values vs Uniform(0,1)"):
    """
    Generates a high-resolution custom QQ plot with specific tick intervals 
    and thin styling for both the data and the reference line.

    :param p_values: Array of p-values (numpy array).
    :param dist: Distribution to compare against (default: Uniform).
    :param title: The title of the plot.
    :return: The matplotlib Figure object.
    """
    
    # 1. Setup the figure with high resolution (DPI=300)
    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)

    # 2. Create the QQ plot using statsmodels
    # We use small markers without edges to make the data look like a thin line
    sm.qqplot(
        p_values, 
        dist=dist, 
        line="45", 
        ax=ax,
        marker='.',         # Simple dot marker
        markersize=2,       # Very small size to mimic a thin line
        markeredgewidth=0,  # Remove marker edges to prevent thickening
        alpha=0.8           # Slight transparency
    )

    # 3. Adjust the thickness of the reference line (the 45-degree line)
    # Loop through all lines in the axes to force a thin linewidth
    for line in ax.get_lines():
        line.set_linewidth(1.0) 

    # 4. Set ticks every 0.1 (from 0 to 1.0)
    ticks = np.arange(0, 1.1, 0.1)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)

    # 5. Titles and Labels
    plt.title(title, fontsize=14)
    plt.xlabel("Expected p-values", fontsize=12)
    plt.ylabel("Observed p-values", fontsize=12)
    
    # Add a subtle grid
    plt.grid(True, linestyle='--', alpha=0.4) 
    
    plt.show()
    
# --- Usage Example ---

# Assuming you have your dataframe:
# p = df_extreme["p_value"].values

# Call the function:
# plot_custom_qq(p)




def plot_log_qq(p_values, title="QQ Plot (-log10)"):
    # 1. Remove NaNs and zeros to avoid log errors
    p_values = p_values[~np.isnan(p_values)]
    p_values = p_values[p_values > 0] 

    # 2. Sort the observed p-values
    p_observed = np.sort(p_values)
    
    # 3. Calculate Expected p-values (Uniform distribution)
    #    Expected under null hypothesis: uniform distribution from 0 to 1
    n = len(p_observed)
    p_expected = np.arange(1, n + 1) / (n + 1)

    # 4. Transform both to -log10
    log_observed = -np.log10(p_observed)
    log_expected = -np.log10(p_expected)

    # 5. Setup the figure
    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)

    # 6. Create the scatter plot
    #    We use 'plot' with linestyle='none' to mimic the previous formatting exactly
    ax.plot(
        log_expected, 
        log_observed, 
        marker='.', 
        linestyle='none', 
        markersize=2, 
        markeredgewidth=0, 
        alpha=0.8
    )

    # 7. Add the 45-degree reference line
    #    It should go from 0 to the maximum value in the plot
    max_val = max(log_expected.max(), log_observed.max())
    ax.plot([0, max_val], [0, max_val], color='red', linewidth=1.0)

    # 8. Titles and Labels (Updated for log scale)
    plt.title(title, fontsize=14)
    plt.xlabel(r"Expected $-\log_{10}(p)$", fontsize=12)
    plt.ylabel(r"Observed $-\log_{10}(p)$", fontsize=12)

    # 9. Dynamic ticks (Standard 0-1 ticks don't work for log scale)
    #    Let pyplot handle it, or set strictly if needed based on max_val
    #    Here we let it be automatic but ensure grid is nice
    plt.grid(True, linestyle='--', alpha=0.4)

    plt.show()

