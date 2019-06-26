"""
Plottint functions shared between the examples.
"""


def remove_ticks(ax):
    ax.tick_params(
        axis='both',          # changes apply to the x-axis
        which='both',       # both major and minor ticks are affected
        bottom=False,       # ticks along the bottom edge are off
        top=False,          # ticks along the top edge are off
        left=False,         # ticks along left edge are off
        labelbottom=False,  # labels along the bottom edge are off
        labelleft=False,    # labels along the left edge are off
        labelright=False)   # labels along the right edge are off
