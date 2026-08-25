"""Shared figure style.

Categorical hues are assigned to ARMS in fixed order and never cycled — the
colour follows the method, not its rank, so a figure that drops an arm does not
repaint the survivors. Palette validated for colour-vision deficiency
(worst adjacent pair ΔE 9.1 protan / 22.9 normal).
"""
import matplotlib as mpl
import matplotlib.pyplot as plt

# fixed categorical slots — do not reorder
ARM_COLOR = {
    'A':         '#2a78d6',   # slot 1  sk2bGrow
    'B':         '#eb6834',   # slot 2  sk2bGrow, Pilea-parity estimator
    'C_default': '#1baf7a',   # slot 3  Pilea, defaults
    'C_relaxed': '#eda100',   # slot 4  Pilea, gates off
    'E':         '#1baf7a',   # FracMinHash sketch + our estimator -- shares the
                              # green of the other FracMinHash arms on purpose:
                              # in the 2x2 the colour encodes the *sketch*, and
                              # the panel encodes the estimator.
}
ARM_LABEL = {
    'A':         'sk2bGrow',
    'B':         'sk2bGrow (Pilea-parity estimator)',
    'C_default': 'Pilea (defaults)',
    'C_relaxed': 'Pilea (gates off)',
    'E':         'FracMinHash sketch + coordinate fit',
}
ARM_ORDER = ['A', 'B', 'C_default', 'C_relaxed', 'E']

#: The attribution 2x2: (estimator, sketch) -> arm id.
ATTRIBUTION = {
    ('coordinate V-fit', '2bRAD anchors'): 'A',
    ('coordinate V-fit', 'FracMinHash'): 'E',
    ('rank regression', '2bRAD anchors'): 'B',
    ('rank regression', 'FracMinHash'): 'C_relaxed',
}
SKETCH_COLOR = {'2bRAD anchors': '#2a78d6', 'FracMinHash': '#1baf7a'}

#: Reference conditions in the fragmentation experiment. A different categorical
#: dimension from ARMS, so it takes the same four validated slots in its own
#: fixed order. `complete` keeps sk2bGrow's blue because it *is* sk2bGrow, on the
#: reference it was designed for.
REF_COLOR = {
    'complete':  '#2a78d6',
    'frag':      '#eb6834',
    'scafRel':   '#1baf7a',
    'pileaFrag': '#eda100',
}
REF_LABEL = {
    'complete':  'complete chromosome',
    'frag':      '100 contigs',
    'scafSelf':  'scaffolded vs itself',
    'scafRel':   'scaffolded vs a relative',
    'pileaFrag': 'Pilea on 100 contigs',
}
REF_ORDER = ['complete', 'frag', 'scafRel', 'pileaFrag']

#: Estimators compared on a fragmented reference (Fig 8d). The V-fit keeps the
#: `frag` orange it already has in a/b/c and Pilea keeps its amber, so only the
#: order-free prototype needs a slot; it takes the unused blue.
EST_COLOR = {'frag': '#eb6834', 'spread_frag': '#2a78d6', 'pileaFrag': '#eda100'}
EST_LABEL = {'frag': 'coordinate V-fit',
             'spread_frag': 'order-free spread MLE',
             'pileaFrag': 'Pilea (sorted rank)'}
EST_ORDER = ['frag', 'spread_frag', 'pileaFrag']

INK, INK2, MUTED = '#1a1a1a', '#4a4a4a', '#8a8a8a'
GRID, SURFACE = '#e6e6e6', '#fcfcfb'


def apply():
    mpl.rcParams.update({
        'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE,
        'savefig.facecolor': SURFACE, 'savefig.bbox': 'tight', 'savefig.dpi': 200,
        'font.family': 'sans-serif',
        'font.sans-serif': ['Helvetica Neue', 'Helvetica', 'Arial', 'DejaVu Sans'],
        'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9,
        'axes.edgecolor': MUTED, 'axes.linewidth': 0.8, 'axes.labelcolor': INK,
        'axes.spines.top': False, 'axes.spines.right': False,   # recessive frame
        'xtick.color': INK2, 'ytick.color': INK2,
        'xtick.labelsize': 8, 'ytick.labelsize': 8,
        'grid.color': GRID, 'grid.linewidth': 0.6,
        'legend.frameon': False, 'legend.fontsize': 8,
        'lines.linewidth': 2.0, 'lines.markersize': 5,           # thin marks
    })


def grid(ax, axis='y'):
    ax.grid(True, axis=axis, zorder=0)
    ax.set_axisbelow(True)
