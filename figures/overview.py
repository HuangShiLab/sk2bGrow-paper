"""Figure 1 — method overview.

Panel (a) is a schematic; (b)-(d) are real output from one benchmark sample
(E. coli K-12 in medium M1, subsampled to 2x), so the figure explains the method
*and* shows it working in the low-coverage regime the paper is about. Everything
is read from `data/exemplar_*`; nothing is recomputed from reads.
"""
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

sys.path.insert(0, str(Path(__file__).parent))
from style import ARM_COLOR, INK, INK2, MUTED, grid

BLUE, ORANGE = ARM_COLOR['A'], ARM_COLOR['B']
GLEN = 4_641_652          # NC_000913.3
EX = 'M1_2x'
TOP = 'CjeI'              # densest enzyme; the one panel b/c follow


# --------------------------------------------------------------------------- a
def _box(ax, x, y, w, h, text, fc, ec, fs=7.2, lw=0.9):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.006,rounding_size=0.012',
                                fc=fc, ec=ec, lw=lw, zorder=2))
    ax.text(x + w / 2, y + h / 2, text, ha='center', va='center', fontsize=fs,
            color=INK, zorder=3, linespacing=1.35)


def _arrow(ax, x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle='-|>', mutation_scale=8,
                                 color=MUTED, lw=1.1, shrinkA=0, shrinkB=0, zorder=1))


def panel_a(ax):
    """Pipeline. Blue-filled boxes are the steps that differ from Pilea."""
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off')
    NEW, PLAIN = '#dce9f9', '#f2f2f0'
    H, Y = 0.42, 0.20

    stages = [
        (0.005, 0.150, 'Reference genomes\n(FASTA / MAG)', PLAIN, False),
        (0.170, 0.185, 'in-silico digest\n16 Type IIB enzymes\n43,735 anchors', NEW, True),
        (0.375, 0.160, 'count reads\ntag scan, ≤2 mismatch\nEM reassignment', PLAIN, False),
        (0.555, 0.185, 'per-enzyme windows\n100 anchors each\nZTP / ZTNB rate ± SE', NEW, True),
        (0.760, 0.235, 'V-shape fit on coordinates\ninverse-variance fusion\nCochran Q, log₂PTR ± CI', NEW, True),
    ]
    for x, w, t, fc, new in stages:
        _box(ax, x, Y, w, H, t, fc, BLUE if new else MUTED, lw=1.1 if new else 0.8)
    for i in range(len(stages) - 1):
        _arrow(ax, stages[i][0] + stages[i][1], Y + H / 2, stages[i + 1][0], Y + H / 2)

    _box(ax, 0.375, 0.80, 0.160, 0.17, 'shotgun reads  /  2bRAD reads', '#ffffff', MUTED, fs=7.0)
    _arrow(ax, 0.455, 0.80, 0.455, Y + H)

    ax.text(0.005, 0.03, 'blue: steps a FracMinHash pipeline cannot take — anchors are motif-defined, so their '
                         'coordinates are known before any read is seen',
            fontsize=6.8, color=INK2, style='italic')
    ax.text(-0.01, 1.02, 'a', fontsize=11, fontweight='bold', color=INK, va='top')


# ----------------------------------------------------------------------- b/c/d
def circ_dist(x, ori):
    d = np.abs(x - ori)
    return np.minimum(d, GLEN - d)


def _v_model(d, p):
    """Plain V in log2 space: log2 mu(d) = a - p * d / (L/2).  This is the model
    fit.py fits; `a` is a nuisance intercept, recovered from the sample mean, so
    the drawn line is the fitted line and not a refit."""
    return -p * d / (GLEN / 2)


def _clip(ax, ys, pad=0.35):
    """Rate windows occasionally drop out entirely (a repeat-masked or deleted
    region), which would otherwise set a 12-unit y-range around a 2-unit signal.
    Clip to the bulk and mark whatever falls off the bottom."""
    lo, hi = np.percentile(ys, [2, 100])
    lo = min(lo, np.percentile(ys, 5) - 0.1)
    ax.set_ylim(lo - pad, hi + pad)
    return lo - pad


def _mark_offscale(ax, xs, ys, floor, color, at='tl'):
    off = ys < floor
    if not off.any():
        return
    ax.scatter(np.asarray(xs)[off], np.full(off.sum(), floor + 0.06), s=22, marker='v',
               color=color, lw=0, zorder=6, clip_on=False)
    x, y, ha, va = (0.02, 0.975, 'left', 'top') if at == 'tl' else (0.98, 0.03, 'right', 'bottom')
    ax.text(x, y, f'{int(off.sum())} window off scale', transform=ax.transAxes,
            fontsize=6.2, color=MUTED, ha=ha, va=va)


def panel_b(ax, w, pe):
    """Window rate against genomic coordinate — what deterministic anchors buy."""
    grid(ax, axis='both')
    ori = float(pe['ori'].iloc[0])
    s = w[w.enzyme == TOP]
    floor = _clip(ax, s.log2_rate.values)

    ax.scatter(w.global_mid / 1e6, w.log2_rate, s=4, c=MUTED, alpha=0.30, lw=0, zorder=2,
               label='all 16 enzymes')
    ax.scatter(s.global_mid / 1e6, s.log2_rate, s=11, c=BLUE, lw=0, zorder=4,
               label=f'{TOP}  (n = {len(s)})')
    _mark_offscale(ax, s.global_mid / 1e6, s.log2_rate.values, floor, BLUE)

    p = float(pe.loc[pe.enzyme == TOP, 'log2_ptr'].iloc[0])
    a = float(s.log2_rate.median()) - float(np.median(_v_model(circ_dist(s.global_mid.values, ori), p)))
    xs = np.linspace(0, GLEN, 800)
    ax.plot(xs / 1e6, _v_model(circ_dist(xs, ori), p) + a, color=BLUE, lw=2.0, zorder=5)
    ax.axvline(ori / 1e6, color=INK2, lw=1.0, ls=(0, (4, 3)), zorder=3)
    ax.text(ori / 1e6 - 0.08, ax.get_ylim()[1], 'ori ', fontsize=7, color=INK2, ha='right', va='top')
    ax.set_xlabel('genome coordinate (Mb)'); ax.set_ylabel('log₂ window rate')
    ax.legend(loc='lower right', handletextpad=0.35, borderpad=0.15, labelspacing=0.25,
              markerscale=2.2)
    ax.set_title('coordinates retained: direct V fit', fontsize=8.6, color=INK, pad=4)
    ax.text(-0.24, 1.14, 'b', transform=ax.transAxes, fontsize=11, fontweight='bold', color=INK, va='top')


def panel_c(ax, w):
    """The same windows with the coordinate thrown away — the rank-regression
    skeleton iRep and Pilea are restricted to."""
    grid(ax, axis='both')
    s = w[w.enzyme == TOP].sort_values('log2_rate')
    n = len(s); r = np.arange(n)
    y = s.log2_rate.values
    floor = _clip(ax, y)

    ax.scatter(r, y, s=11, c=ORANGE, lw=0, zorder=4)
    b, a = np.polyfit(r, y, 1)                     # slope x n windows = log2 PTR
    ax.plot(r, a + b * r, color=ORANGE, lw=2.0, zorder=5, label=f'all windows: {b * n:.2f}')
    _mark_offscale(ax, r, y, floor, ORANGE, at='br')

    keep = y > floor                               # drop the single dropout window
    b2, a2 = np.polyfit(r[keep], y[keep], 1)
    ax.plot(r[keep], a2 + b2 * r[keep], color=ORANGE, lw=1.4, ls=(0, (4, 2.5)), zorder=5,
            label=f'minus 1 window: {b2 * keep.sum():.2f}')
    ax.legend(loc='upper left', title='log₂(PTR) =', title_fontsize=6.8, fontsize=6.8,
              handletextpad=0.5, labelspacing=0.25, borderpad=0.2, alignment='left')
    ax.set_xlabel('window rank'); ax.set_ylabel('log₂ window rate')
    ax.set_title('coordinates discarded: rank regression', fontsize=8.6, color=INK, pad=4)
    ax.text(-0.24, 1.14, 'c', transform=ax.transAxes, fontsize=11, fontweight='bold', color=INK, va='top')


def panel_d(ax, pe, out):
    """Forest plot: 16 independent measurement channels, then fusion."""
    pe = pe.sort_values('n_anchors').reset_index(drop=True)
    y = np.arange(len(pe))
    fused = float(out['log2(PTR)'].iloc[0])
    lo, hi = float(out['ci_low'].iloc[0]), float(out['ci_high'].iloc[0])

    ax.axvspan(lo, hi, color=BLUE, alpha=0.18, lw=0, zorder=1)
    ax.axvline(fused, color=BLUE, lw=1.8, zorder=2)
    ax.errorbar(pe.log2_ptr, y, xerr=1.96 * pe.se, fmt='o', ms=3.6, lw=0, elinewidth=1.1,
                color=INK2, ecolor=MUTED, zorder=4)
    ax.set_yticks(y); ax.set_yticklabels(pe.enzyme, fontsize=6.4)
    ax.set_ylim(-0.8, len(pe) - 0.2)
    ax.set_xlabel('log₂(PTR) per enzyme')
    ax.tick_params(axis='y', length=0)
    ax.grid(True, axis='x', zorder=0); ax.set_axisbelow(True)
    ax.text(0.0, -0.235, f'fused {fused:.2f}  [{lo:.2f}, {hi:.2f}]', transform=ax.transAxes,
            fontsize=7.0, color=BLUE, ha='left', va='center')
    ax.text(0.0, -0.325, f'I² = {float(out["enzyme_i2"].iloc[0]):.2f}, random effects',
            transform=ax.transAxes, fontsize=6.6, color=INK2, ha='left', va='center')
    ax.set_title('16 enzymes = 16 strata', fontsize=8.6, color=INK, pad=4)
    ax.text(-0.46, 1.14, 'd', transform=ax.transAxes, fontsize=11, fontweight='bold', color=INK, va='top')


def build(root):
    root = Path(root)
    w = pd.read_csv(root / 'data' / f'exemplar_{EX}_windows.rates.tsv', sep='\t')
    pe = pd.read_csv(root / 'data' / f'exemplar_{EX}_per_enzyme.tsv', sep='\t')
    out = pd.read_csv(root / 'data' / f'exemplar_{EX}_output.tsv', sep='\t')

    fig = plt.figure(figsize=(7.2, 4.8))
    gs = fig.add_gridspec(2, 3, height_ratios=[0.60, 1.0], hspace=0.46, wspace=0.42,
                          left=0.085, right=0.985, top=0.975, bottom=0.135)
    panel_a(fig.add_subplot(gs[0, :]))
    panel_b(fig.add_subplot(gs[1, 0]), w, pe)
    panel_c(fig.add_subplot(gs[1, 1]), w)
    panel_d(fig.add_subplot(gs[1, 2]), pe, out)
    return fig
