# Shared plotting constants and helpers for the analysis notebooks.
#
# This is a trimmed copy of the lab's personal `const.py`: it keeps only the
# names actually used by the notebooks in this repository (figure styling,
# colors, and the save_fig / set_plot_style / set_equal_plot_limits helpers).
# The original module's large `MPRA_data_paths` dictionary was removed because
# no notebook here uses it.

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np


# Standard font size and family
FONT_SIZE_small = 18
FONT_SIZE_big = 20

FONT_FAMILY = 'arial'

# Colors for controls / highlights
pos_active_ctrl_color = 'g'
neg_active_ctrl_color = 'r'
highlight_color = 'y'

DIFF_ACTIVITY_COLOR = 'green'

# Custom colormaps for scatterplots
colors = ["#EBF4FF", "#E1ECFA", "#D0E0F5", "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026"]
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", colors, N=256)

colors = ["#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026"]
custom_cmap_bolder = LinearSegmentedColormap.from_list("custom_cmap", colors, N=256)


def set_plot_style():
    """
    Set standardized figure settings for matplotlib.
    """
    plt.rcParams.update({
        'axes.titlesize': FONT_SIZE_big,
        'axes.labelsize': FONT_SIZE_big,
        'xtick.labelsize': FONT_SIZE_small,
        'ytick.labelsize': FONT_SIZE_small,
        'legend.fontsize': FONT_SIZE_small,
        'legend.title_fontsize': FONT_SIZE_small,
        'axes.labelweight': 'bold',
#        'font.family': FONT_FAMILY, # TODO: change font
        'axes.linewidth': 1.0,
        'figure.figsize': (8, 8),
        'axes.grid': False,  # No grid
        'axes.spines.top': False,  # Top border off
        'axes.spines.right': False,  # Right border off
        'figure.facecolor': 'none',
        'axes.facecolor': 'none',
#        'legend.frameon': False  # No frame for legend
    })


def save_fig(fig, name, path):
    """
    Save the figure to the specified path in PNG and SVG formats.
    """
    fig.savefig(f"{path}/{name}.png", dpi=500, bbox_inches='tight', transparent=True)
    fig.savefig(f"{path}/{name}.svg", bbox_inches="tight", transparent=True)


def set_equal_plot_limits(x, y):
    """
    Sets the x and y axis limits to the same range based on the min and max values of x and y.

    Parameters:
    x (array-like): Data for the x-axis.
    y (array-like): Data for the y-axis.
    """
    min_limit = min(np.min(x), np.min(y))
    max_limit = max(np.max(x), np.max(y))

    plt.xlim([min_limit, max_limit])
    plt.ylim([min_limit, max_limit])


plot_color_pallete = {
    "default_color": "#AEAEAE",
    "cCRE": "#3D9F95",
    "barcode": "#227C9D",
    "read": "#FFC25F",
    "cCRE-barcode-pair": "#9383B8",
}
