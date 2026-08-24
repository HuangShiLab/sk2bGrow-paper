#!/usr/bin/env python3
"""Regenerate every figure from the committed benchmark tables.

    python3 figures/make_figures.py

Reads only `data/*.tsv`; writes `figures/out/*.{png,pdf}`. No network, no
recomputation — the figures are a pure function of the committed data, so a
reviewer can reproduce them without rerunning the benchmark.
"""
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from style import ARM_COLOR, ARM_LABEL, ARM_ORDER, INK, INK2, MUTED, GRID as GRIDC, apply, grid

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'figures' / 'out'
OUT.mkdir(parents=True, exist_ok=True)
apply()

res = pd.read_csv(ROOT / 'data' / 'results_raw.tsv', sep='\t')
grow = res[res['medium'] != 'RUN_OUT'].copy()
ctl = res[res['medium'] == 'RUN_OUT'].copy()


def save(fig, name):
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'{name}.{ext}')
    plt.close(fig)
    print(f'  wrote {name}.png / .pdf')


# --- Fig 1: accuracy against measured growth rate, across coverage -----------
def fig1():
    """Pilea's two arms are identical wherever both are defined, so they are one
    series here; marker fill encodes whether its shipped gates would report."""
    fig, ax = plt.subplots(figsize=(5.4, 3.4))
    grid(ax)
    covs = sorted(grow['cov'].unique())

    def curve(arm):
        xs, ys = [], []
        for c in covs:
            s = grow[(grow['cov'] == c) & (grow['arm'] == arm) & np.isfinite(grow['log2ptr'])]
            # A constant output (Pilea returns PTR=1.0 for every sample at 0.5x)
            # has no defined correlation. Omit rather than plot a zero, which
            # would read as "scored badly" instead of "no signal at all".
            if len(s) < 3 or s['log2ptr'].nunique() < 2:
                continue
            xs.append(c); ys.append(stats.pearsonr(s['growth_rate'], s['log2ptr'])[0])
        return xs, ys

    for arm in ['A', 'B']:
        xs, ys = curve(arm)
        ax.plot(xs, ys, 'o-', color=ARM_COLOR[arm], label=ARM_LABEL[arm],
                markeredgecolor='white', markeredgewidth=0.8, zorder=3)

    xs, ys = curve('C_relaxed')
    reports = set(grow[(grow['arm'] == 'C_default')].dropna(subset=['log2ptr'])['cov'])
    ax.plot(xs, ys, '-', color=ARM_COLOR['C_relaxed'], label='Pilea', zorder=3)
    for x, y in zip(xs, ys):
        filled = x in reports
        ax.plot([x], [y], 'o', color=ARM_COLOR['C_relaxed'] if filled else 'white',
                markeredgecolor=ARM_COLOR['C_relaxed'], markeredgewidth=1.6, zorder=4)
    ax.plot([], [], 'o', color='white', markeredgecolor=ARM_COLOR['C_relaxed'],
            markeredgewidth=1.6, label='  (open: its own gates refuse)')

    ax.axhline(0.9764, ls=':', lw=1.2, color=MUTED, zorder=1)
    ax.text(10.6, 0.9764, "Pilea, published\n(full depth)", fontsize=7,
            color=MUTED, va='center', ha='left')
    ax.set_xscale('log')
    ax.set_xticks(covs); ax.set_xticklabels([f'{c:g}×' for c in covs])
    ax.set_xlabel('subsampled coverage')
    ax.set_ylabel('Pearson r vs measured growth rate')
    ax.set_ylim(0.3, 1.02); ax.set_xlim(0.42, 13)
    ax.set_title('Accuracy against an independent ground truth', loc='left', color=INK)
    ax.legend(loc='lower right')
    fig.text(0.0, -0.05,
             'Zheng et al. 2020 E. coli K-12, 16 growth media (λ = 0.40–1.72 h⁻¹). '
             'Pilea at its shipped defaults\nreturns no estimate below 10×; at 0.5× with gates off it '
             'returns PTR = 1.0 for every sample, so\nneither has a defined correlation there.',
             fontsize=7, color=MUTED, va='top')
    save(fig, 'fig1_accuracy_vs_coverage')


# --- Fig 2: is the magnitude right, not just the ranking? --------------------
def fig2():
    covs = [1.0, 10.0]
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.5), sharex=True, sharey=True)
    for ax, c in zip(axes, covs):
        grid(ax, axis='both')
        lo, hi = 0.3, 2.1
        ax.plot([lo, hi], [lo, hi], ls='--', lw=1, color=MUTED, zorder=1)
        for arm in ['A', 'C_relaxed']:
            s = grow[(grow['cov'] == c) & (grow['arm'] == arm)].dropna(subset=['log2ptr'])
            if s.empty or s['log2ptr'].nunique() < 2:
                continue
            ax.scatter(s['pred_log2ptr'], s['log2ptr'], s=34, alpha=0.85,
                       color=ARM_COLOR[arm], label=ARM_LABEL[arm],
                       edgecolor='white', linewidth=0.8, zorder=3)
            sl, ic = np.polyfit(s['pred_log2ptr'], s['log2ptr'], 1)
            xx = np.linspace(lo, hi, 10)
            ax.plot(xx, sl * xx + ic, lw=1.4, color=ARM_COLOR[arm], alpha=0.7, zorder=2)
            ax.text(0.05, 0.93 if arm == 'A' else 0.84, f'slope {sl:.2f}',
                    transform=ax.transAxes, fontsize=7.5, color=ARM_COLOR[arm])
        ax.set_title(f'{c:g}× coverage', loc='left', color=INK)
        ax.set_xlabel('predicted log₂(PTR)   (Zheng λC/ln2)')
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect('equal')
    axes[0].set_ylabel('estimated log₂(PTR)')
    axes[0].legend(loc='lower right')
    fig.text(0.0, -0.04,
             'Dashed line is y = x. A slope below 1 means the PTR range is compressed even where '
             'the ranking\nis correct — which a correlation coefficient hides. The predicted value '
             'is not independent: Zheng\nderived C from this same sequencing by marker-frequency analysis.',
             fontsize=7, color=MUTED, va='top')
    save(fig, 'fig2_magnitude')


# --- Fig 3: the negative control --------------------------------------------
def fig3():
    """Pilea's two arms coincide wherever both are defined, so they are drawn as
    one series (as in fig 1) rather than as a phantom extra legend key."""
    from matplotlib.lines import Line2D
    fig, ax = plt.subplots(figsize=(5.6, 3.3))
    grid(ax)
    covs = sorted(ctl['cov'].unique())
    keys = []
    for arm, lab in [('A', ARM_LABEL['A']), ('B', ARM_LABEL['B']), ('C_relaxed', 'Pilea')]:
        s_ = ctl[ctl['arm'] == arm].dropna(subset=['log2ptr']).sort_values('cov')
        if s_.empty:
            continue
        ax.plot(s_['cov'], s_['log2ptr'], 'o-', color=ARM_COLOR[arm],
                markeredgecolor='white', markeredgewidth=0.8, zorder=3)
        keys.append(Line2D([], [], marker='o', color=ARM_COLOR[arm], ms=6,
                           markeredgecolor='white', markeredgewidth=0.8, label=lab))
    ax.axhline(0, lw=1.2, color=MUTED, zorder=1)
    ax.annotate('truth: not growing', xy=(2.0, 0), xytext=(2.0, -0.34),
                fontsize=7.5, color=MUTED, ha='center',
                arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.8))
    ax.set_xscale('log')
    ax.set_xticks(covs); ax.set_xticklabels([f'{c:g}×' for c in covs])
    ax.set_ylim(-0.55, 2.45)
    ax.set_xlabel('subsampled coverage')
    ax.set_ylabel('estimated log₂(PTR)')
    ax.set_title('Negative control: a stationary-phase culture', loc='left', color=INK)
    ax.legend(handles=keys, loc='upper right')
    fig.text(0.0, -0.06,
             'The sorted-regression estimator reports log₂(PTR) = 2.17 at 0.5× for a culture that is not '
             'growing —\nit manufactures a gradient out of rank-ordered noise. The coordinate fit stays '
             'within 0.11 of zero at\nevery coverage. Pilea reports only at 10× under its own gates; the '
             'curve shown has them disabled.',
             fontsize=7, color=MUTED, va='top')
    save(fig, 'fig3_negative_control')


# --- Fig 4: what is responsible — sketch or estimator? ----------------------
def fig4():
    """Dot plot, not bars: every value sits in 0.45-0.98, so a zero-based bar
    chart wastes half the panel. A dot plot legitimately takes a truncated axis
    where a bar chart may not."""
    from matplotlib.lines import Line2D
    fig, ax = plt.subplots(figsize=(6.0, 3.4))
    grid(ax, axis='x')
    covs = [0.5, 1.0, 2.0, 5.0]
    series = [('A', 'anchors + coordinate fit'),
              ('B', 'anchors + sorted regression'),
              ('C_relaxed', 'FracMinHash sketch + Pilea')]
    for row, c in enumerate(covs):
        vals = {}
        for arm, _ in series:
            s_ = grow[(grow['cov'] == c) & (grow['arm'] == arm)].dropna(subset=['log2ptr'])
            if len(s_) >= 3 and s_['log2ptr'].nunique() > 1:
                vals[arm] = stats.pearsonr(s_['growth_rate'], s_['log2ptr'])[0]
        if len(vals) > 1:
            ax.plot([min(vals.values()), max(vals.values())], [row, row],
                    color=GRIDC, lw=3, solid_capstyle='round', zorder=2)
        # Alternate label side when two dots crowd, so numbers never overlap.
        placed = []
        for arm, _ in series:
            if arm not in vals:
                continue
            v = vals[arm]
            ax.plot([v], [row], 'o', ms=9, color=ARM_COLOR[arm],
                    markeredgecolor='white', markeredgewidth=1.2, zorder=3)
            below = any(abs(v - q) < 0.045 for q in placed)
            ax.text(v, row + (0.30 if below else -0.22), f'{v:.2f}', ha='center',
                    va='top' if below else 'bottom', fontsize=7, color=INK2)
            placed.append(v)
    ax.set_yticks(range(len(covs)))
    ax.set_yticklabels([f'{c:g}×' for c in covs])
    ax.set_xlim(0.38, 1.04)
    ax.set_ylim(len(covs) - 0.5, -0.5)
    ax.set_xlabel('Pearson r vs measured growth rate')
    ax.set_ylabel('subsampled coverage')
    ax.set_title('The estimator, not the sketch, carries the gain', loc='left', color=INK)
    # explicit proxies: a series absent from the first row still needs a key
    ax.legend(handles=[Line2D([], [], marker='o', ls='', ms=8, color=ARM_COLOR[a],
                              markeredgecolor='white', markeredgewidth=1.2, label=l)
                       for a, l in series],
              loc='upper center', bbox_to_anchor=(0.5, -0.28), ncol=3,
              columnspacing=1.2, handletextpad=0.4)
    fig.text(0.0, -0.30,
             'Holding the sketch fixed and swapping the estimator (the two anchor rows) moves accuracy '
             'far more\nthan holding the estimator fixed and swapping the sketch (sorted regression, '
             'anchors vs FracMinHash).\nUnder a sorted-regression estimator the deterministic anchors '
             'are behind FracMinHash at 1× — 0.61 vs 0.89.',
             fontsize=7, color=MUTED, va='top')
    save(fig, 'fig4_attribution')


if __name__ == '__main__':
    print('regenerating figures ->', OUT)
    fig1(); fig2(); fig3(); fig4()
    print('done')
