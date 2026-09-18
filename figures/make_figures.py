#!/usr/bin/env python3
"""Regenerate every figure from the committed benchmark tables.

    python3 figures/make_figures.py

Reads only `data/*.tsv`; writes `figures/out/*.{png,pdf}`. No network, no
recomputation — the figures are a pure function of the committed data, so a
reviewer can reproduce them without rerunning the benchmark.
"""
from pathlib import Path
import os
import sys

# Matplotlib stamps the wall clock into every PDF, so regenerating dirties all
# of them even when nothing changed. Pin it: the figures are a pure function of
# the committed data and their bytes should be too.
os.environ.setdefault('SOURCE_DATE_EPOCH', '0')

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patheffects as pe

sys.path.insert(0, str(Path(__file__).parent))
from style import (ARM_COLOR, ARM_LABEL, ARM_ORDER, ATTRIBUTION, INK, INK2, MUTED,
                   SKETCH_COLOR, REF_COLOR, REF_LABEL, REF_ORDER,
                   EST_COLOR, EST_LABEL, EST_ORDER,
                   GRID as GRIDC, SURFACE, apply, grid)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'figures' / 'out'
OUT.mkdir(parents=True, exist_ok=True)
apply()

res = pd.read_csv(ROOT / 'data' / 'results_raw.tsv', sep='\t')
grow = res[res['medium'] != 'RUN_OUT'].copy()
ctl = res[res['medium'] == 'RUN_OUT'].copy()


NL = chr(10)


def save(fig, name):
    for ext in ('png', 'pdf'):
        fig.savefig(OUT / f'{name}.{ext}')
    plt.close(fig)
    print(f'  wrote {name}.png / .pdf')


# --- Fig 1: method overview -------------------------------------------------
def fig1():
    """Schematic plus three real panels from one 2x sample; see overview.py."""
    import overview
    save(overview.build(ROOT), 'fig1_overview')


#: Fig 2's caption: the three standalone captions it merges, panel-prefixed.
#: The middle paragraph describes both scatter panels, so it carries both
#: letters.
CAP2_HEAD = (
    'a, Zheng et al. 2020 E. coli K-12, 16 growth media (\u03bb = 0.40\u20131.72 h\u207b\u00b9). '
    'Pilea at its shipped defaults\nreturns no estimate below 10\u00d7; at 0.5\u00d7 with '
    'gates off it returns PTR = 1.0 for every sample, so\nneither has a defined '
    'correlation there.'
    + NL + NL +
    'b, c, Dashed line is y = x. A slope below 1 means the PTR range is compressed '
    'even where the ranking\nis correct \u2014 which a correlation coefficient hides. '
    'The predicted value is not independent: Zheng\nderived C from this same '
    'sequencing by marker-frequency analysis.'
    + NL + NL +
    'd, The sorted-regression estimator reports log\u2082(PTR) = {b_ro:.2f} at 0.5\u00d7 for '
    'a culture that is not growing \u2014\nit manufactures a gradient out of '
    'rank-ordered noise. The coordinate fit stays within {a_ro:.2f} of zero at\n'
    'every coverage. Pilea reports only at 10\u00d7 under its own gates; the curve '
    'shown has them disabled.'
)


def cap2():
    """Panel-d numbers computed from the RUN_OUT control so the caption cannot
    drift from the data (it quotes two specific values)."""
    b_ro = ctl[(ctl['arm'] == 'B') & (ctl['cov'] == sorted(ctl['cov'].unique())[0])]['log2ptr'].iloc[0]
    a_ro = ctl[ctl['arm'] == 'A']['log2ptr'].abs().max()
    return CAP2_HEAD.format(b_ro=b_ro, a_ro=a_ro)


# --- Fig 2: Zheng benchmark — accuracy, magnitude, negative control ----------
def fig2():
    """The E. coli K-12 benchmark as one row: accuracy against coverage (a),
    estimated vs predicted magnitude at 1x and 10x (b, c), and the stationary
    negative control (d)."""
    covs = sorted(grow['cov'].unique())
    fig = plt.figure(figsize=(12.4, 3.45))
    gs = fig.add_gridspec(1, 3, width_ratios=[2.00, 2.208 * 2 + 0.10, 2.55],
                          wspace=0.48, left=0.03, right=0.985,
                          top=0.795, bottom=0.155)
    gs_m = gs[0, 1].subgridspec(1, 2, wspace=0.045, width_ratios=[1, 1])
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs_m[0, 0])
    ax_c = fig.add_subplot(gs_m[0, 1], sharex=ax_b, sharey=ax_b)
    ax_d = fig.add_subplot(gs[0, 2])

    # -- a: accuracy against measured growth rate, across coverage (old fig 2)
    grid(ax_a)

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
        ax_a.plot(xs, ys, 'o-', color=ARM_COLOR[arm], label=ARM_LABEL[arm],
                  markeredgecolor='white', markeredgewidth=0.8, zorder=3)

    xs, ys = curve('C_relaxed')
    reports = set(grow[(grow['arm'] == 'C_default')].dropna(subset=['log2ptr'])['cov'])
    ax_a.plot(xs, ys, '-', color=ARM_COLOR['C_relaxed'], label='Pilea', zorder=3)
    for x, y in zip(xs, ys):
        filled = x in reports
        ax_a.plot([x], [y], 'o', color=ARM_COLOR['C_relaxed'] if filled else 'white',
                  markeredgecolor=ARM_COLOR['C_relaxed'], markeredgewidth=1.6, zorder=4)
    ax_a.plot([], [], 'o', color='white', markeredgecolor=ARM_COLOR['C_relaxed'],
              markeredgewidth=1.6, label='  (open: its own gates refuse)')

    ax_a.axhline(0.9764, ls=':', lw=1.2, color=MUTED, zorder=1)
    # The reference line's label goes above the line at the far left: to its
    # right the 5x/10x points crowd the line, and below it arm A runs.
    ax_a.text(0.45, 0.985, "Pilea, published\n(full depth)", fontsize=7,
              color=MUTED, va='bottom', ha='left')
    ax_a.set_xscale('log')
    ax_a.set_xticks(covs); ax_a.set_xticklabels([f'{c:g}×' for c in covs])
    ax_a.xaxis.set_minor_locator(mticker.NullLocator())
    ax_a.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax_a.set_xlabel('subsampled coverage')
    ax_a.set_ylabel('Pearson r vs measured growth rate')
    # arm B collapses to r = 0.16 at 0.5x; the axis has to reach it, or the
    # series would leave the frame with no indication that it did.
    ax_a.set_ylim(0.1, 1.03); ax_a.set_xlim(0.42, 13)
    ax_a.set_title('Accuracy against an independent\nground truth', loc='left',
                   color=INK, fontsize=9, pad=14)
    ax_a.legend(loc='lower right', fontsize=7.5)

    # -- b, c: is the magnitude right, not just the ranking? (old fig 3) ------
    for ax, c in zip([ax_b, ax_c], [1.0, 10.0]):
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
        ax.set_title(f'{c:g}× coverage', loc='left', color=INK, fontsize=9, pad=14)
        ax.set_xlabel('predicted log₂(PTR)   (Zheng λC/ln2)')
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect('equal')
    ax_b.set_ylabel('estimated log₂(PTR)')
    ax_b.legend(loc='lower right', fontsize=7.5)
    plt.setp(ax_c.get_yticklabels(), visible=False)

    # -- d: the negative control (old fig 4) ---------------------------------
    from matplotlib.lines import Line2D
    grid(ax_d)
    d_covs = sorted(ctl['cov'].unique())
    keys = []
    for arm, lab in [('A', ARM_LABEL['A']), ('B', ARM_LABEL['B']), ('C_relaxed', 'Pilea')]:
        s_ = ctl[ctl['arm'] == arm].dropna(subset=['log2ptr']).sort_values('cov')
        if s_.empty:
            continue
        ax_d.plot(s_['cov'], s_['log2ptr'], 'o-', color=ARM_COLOR[arm],
                  markeredgecolor='white', markeredgewidth=0.8, zorder=3)
        keys.append(Line2D([], [], marker='o', color=ARM_COLOR[arm], ms=6,
                           markeredgecolor='white', markeredgewidth=0.8, label=lab))
    ax_d.axhline(0, lw=1.2, color=MUTED, zorder=1)
    ax_d.annotate('truth: not growing', xy=(2.0, 0), xytext=(2.0, -0.34),
                  fontsize=7.5, color=MUTED, ha='center',
                  arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.8))
    ax_d.set_xscale('log')
    ax_d.set_xticks(d_covs); ax_d.set_xticklabels([f'{c:g}×' for c in d_covs])
    ax_d.set_ylim(-0.55, 2.45)
    ax_d.set_xlabel('subsampled coverage')
    ax_d.set_ylabel('estimated log₂(PTR)')
    ax_d.set_title('Negative control: a\nstationary-phase culture', loc='left',
                   color=INK, fontsize=9, pad=14)
    ax_d.legend(handles=keys, loc='upper right', fontsize=7.5)

    for ax, letter, dx in ((ax_a, 'a', -0.36), (ax_b, 'b', -0.30),
                           (ax_c, 'c', -0.15), (ax_d, 'd', -0.26)):
        ax.text(dx, 1.10, letter, transform=ax.transAxes, fontsize=11,
                fontweight='bold', color=INK, va='top')
    fig.text(0.0, -0.02, cap2(), fontsize=7, color=MUTED, va='top')
    save(fig, 'fig2_zheng_benchmark')


# --- Fig 3: multi-strain simulation — accuracy, recall and cost -------------
def fig3():
    """Three panels because three different things matter and they trade off:
    how often a method answers, how right it is, and what it costs."""
    f = ROOT / 'data' / 'sim_results.tsv'
    if not f.exists():
        print('  (skipping fig3: sim_results.tsv absent)'); return
    d = pd.read_csv(f, sep='\t')
    SIM = {'sk2bGrow': ARM_COLOR['A'],
           'Pilea (defaults)': ARM_COLOR['C_default'],
           'Pilea (gates off)': ARM_COLOR['C_relaxed']}
    covs = sorted(d['coverage'].unique())
    fig, axes = plt.subplots(1, 3, figsize=(10.4, 3.2))

    for ax, (col, lab, lo) in zip(axes, [
            ('recall', 'strains reported (recall)', 0),
            ('rmse', 'RMSE of log₂(PTR)', 0),
            ('seconds', 'wall-clock seconds', 0)]):
        grid(ax)
        for arm, colr in SIM.items():
            s = d[d['arm'] == arm].groupby('coverage')[col].mean().reindex(covs)
            # A lone finite point draws no line, so mark points explicitly too —
            # otherwise Pilea-at-defaults vanishes from the RMSE panel entirely.
            ax.plot(s.index, s.values, '-', color=colr, zorder=3)
            ok = s.dropna()
            ax.plot(ok.index, ok.values, 'o', color=colr, label=arm,
                    markeredgecolor='white', markeredgewidth=0.8, zorder=4)
        ax.set_xscale('log'); ax.set_xticks(covs)
        ax.set_xticklabels([f'{c:g}×' for c in covs])
        ax.xaxis.set_minor_locator(mticker.NullLocator())   # kill "3 x 10^0"
        ax.xaxis.set_minor_formatter(mticker.NullFormatter())
        ax.set_xlabel('coverage per strain')
        ax.set_ylabel(lab)
        ax.set_ylim(bottom=lo)
    axes[0].set_ylim(-0.05, 1.08)
    # The band between the recall = 1 plateau and Pilea-defaults' late rise
    # (x ~ [1, 2.5], y ~ [0.35, 0.75]) is empty; lower right is not, the green
    # curve climbs straight through it.
    axes[0].legend(loc='center left', bbox_to_anchor=(0.03, 0.52), fontsize=7.5)
    axes[0].set_title('answers at all', loc='left', color=INK)
    axes[1].set_title('when it answers, how right', loc='left', color=INK)
    axes[2].set_title('what it costs', loc='left', color=INK)
    fig.text(0.0, -0.08,
             'Multi-strain communities: 16 reference genomes, 4/8/16 strains per sample, V-shaped '
             'profiles, log₂PTR ~ U[0,2].\nRecall and RMSE must be read together — Pilea at its '
             'shipped defaults earns a flattering RMSE by answering only\n22% of cases, the easiest '
             'ones. sk2bGrow answers every case; cost crosses over at 4×, below which it is the '
             'cheaper of the two\narms that answer.',
             fontsize=7, color=MUTED, va='top')
    save(fig, 'fig3_simulation')


#: Fig 4's caption, kept at module level so the hard line breaks stay visible.
CAP4 = (
    'a, the full 2\u00d72 at the Pilea operating point — scale 250, a landmark density '
    '2.4\u00d7 below the enzyme panel — which is where the sketch penalty lives; the '
    'interaction is a fact about this operating point, not a property of the' + NL +
    'sketch. b, the same comparison after density matching (FracMinHash scale 104 '
    'vs the 16-enzyme panel; 9,645 vs 9,422 landmarks/Mb): at \u22651\u00d7 the two landmark '
    'modes are tied at every depth, and at 0.5\u00d7 the sketch scores higher on the' + NL +
    'cells it survives (r 0.82\u20130.84, n = 10\u201311, vs 0.50\u20130.57, n = 15\u201316) — but it '
    'loses 6 of 17 cells to `no downhill origin`. c, the mechanism at 0.5\u00d7: '
    'windows are populated almost identically in the two' + NL +
    'modes (median 17.5 vs 14.1 landmarks/window; detected fraction 0.175 vs '
    '0.185) — deterministic anchors are not what keeps windows populated. The '
    'panel\u2019s low-coverage robustness comes from redundancy across its 16' + NL +
    'strata: a single-strata sketch dies outright when its one stratum fits no '
    'gradient and no fusion can rescue it, while the panel still answers if any '
    'of its 16 strata do.'
)


# --- Fig 4: what is responsible — sketch or estimator? ----------------------
def fig4():
    """Panel a: the full 2x2 at the Pilea operating point. Panel b: the same
    comparison once sketch density is matched to the panel. Panel c: the
    mechanism — per-window landmark counts, which turn out to be identical
    across modes, so the difference is strata redundancy, not anchor kind.
    """
    f1 = ROOT / 'data' / 'f1_sketch' / 'F1_results.tsv'
    f1m = ROOT / 'data' / 'f1_sketch' / 'F1_mechanism.tsv'
    fig = plt.figure(figsize=(11.8, 3.5))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1.02, 0.92], wspace=0.34)
    ax_a1 = fig.add_subplot(gs[0, 0])
    ax_a2 = fig.add_subplot(gs[0, 1], sharey=ax_a1)
    ax_b = fig.add_subplot(gs[0, 2], sharey=ax_a1)
    ax_c = fig.add_subplot(gs[0, 3])
    covs = sorted(grow['cov'].unique())
    ESTIMATORS = [('coordinate V-fit', 'coordinate V-fit (ours)'),
                  ('rank regression', 'rank regression (iRep / Pilea)')]

    def r_at(arm, c):
        s_ = grow[(grow['cov'] == c) & (grow['arm'] == arm)].dropna(subset=['log2ptr'])
        if len(s_) < 3 or s_['log2ptr'].nunique() < 2:
            return np.nan
        return stats.pearsonr(s_['growth_rate'], s_['log2ptr'])[0]

    for ax, (est, title) in zip([ax_a1, ax_a2], ESTIMATORS):
        grid(ax, axis='both')
        drawn = {}
        for sketch, colour in SKETCH_COLOR.items():
            arm = ATTRIBUTION[(est, sketch)]
            ys = [r_at(arm, c) for c in covs]
            ok = [(c, y) for c, y in zip(covs, ys) if np.isfinite(y)]
            if not ok:
                continue
            drawn[sketch] = dict(ok)
            ax.plot([c for c, _ in ok], [y for _, y in ok], 'o-', color=colour,
                    ms=5, markeredgecolor='white', markeredgewidth=1.0, label=sketch)
        # Contrast of the green against the surface is under 3:1, so each series
        # carries a visible label rather than relying on hue alone. Both curves
        # converge at high coverage, so the labels go where the two are furthest
        # apart -- otherwise they collide exactly where the lines do. Under the
        # current data the 0.5x curves sit ~0.06 r apart, so the labels need
        # more vertical room than the original 9/15pt offsets gave them, and
        # when both curves leave the point the same way their labels would
        # share one side — the first series then names its point from above,
        # clear of the other's rise.
        shared = [c for c in covs if all(c in d for d in drawn.values())]
        if shared and len(drawn) == 2:
            d1, d2 = drawn.values()
            cx = max(shared, key=lambda c: abs(d1[c] - d2[c]))
            first = cx == covs[0]
            nxt = [c for c in covs if c > cx]
            items = list(drawn.items())
            sides = [-1 if (nxt and dd.get(nxt[0], dd[cx]) > dd[cx]) else 1
                     for _, dd in items]
            if sides[0] == sides[1]:
                sides = [1, -1]
            for (sketch, dd), side in zip(items, sides):
                # Sit the label on the side the curve is leaving, so a rising
                # series never has its own line drawn through its name.
                ax.annotate(sketch, (cx, dd[cx]), textcoords='offset points',
                            xytext=(4 if first else 0, -21 if side < 0 else 13),
                            ha='left' if first else 'center',
                            fontsize=7, color=SKETCH_COLOR[sketch],
                            path_effects=[pe.withStroke(linewidth=2.5,
                                                        foreground=SURFACE)])
        ax.set_xscale('log')
        ax.set_xticks(covs); ax.set_xticklabels([f'{c:g}×' for c in covs])
        ax.xaxis.set_minor_locator(mticker.NullLocator())
        ax.xaxis.set_minor_formatter(mticker.NullFormatter())
        ax.set_xlabel('subsampled coverage')
        ax.set_title(title, fontsize=9, color=INK, pad=6)
    ax_a1.set_ylabel('Pearson r vs measured growth rate')
    ax_a1.set_ylim(0.1, 1.04)
    ax_a1.text(0.03, 0.03, 'operating point: scale 250,\n2.4× below panel density',
               transform=ax_a1.transAxes, fontsize=6.6, color=MUTED, va='bottom')

    at1 = {k: r_at(v, 1.0) for k, v in ATTRIBUTION.items()}
    gain = {sk: at1[('coordinate V-fit', sk)] - at1[('rank regression', sk)]
            for sk in SKETCH_COLOR}
    inter = gain['2bRAD anchors'] - gain['FracMinHash']

    # --- b: density-matched --------------------------------------------------
    d = pd.read_csv(f1, sep='\t')
    grid(ax_b, axis='both')
    pair = [('A_k16', 'enzyme panel k16', SKETCH_COLOR['2bRAD anchors']),
            ('E_s104', 'FracMinHash s104', SKETCH_COLOR['FracMinHash'])]
    for arm, lab, colour in pair:
        s = d[d['arm'] == arm].sort_values('cov')
        ax_b.plot(s['cov'], s['r'], 'o-', color=colour, ms=5,
                  markeredgecolor='white', markeredgewidth=1.0, label=lab)
    lo = d[(d['arm'] == 'A_k16') & (d['cov'] == covs[0])].iloc[0]
    hi = d[(d['arm'] == 'E_s104') & (d['cov'] == covs[0])].iloc[0]
    # At 0.5x the sketch r is computed on the cells that survived; the dropped
    # cells are the whole story, so both n's are printed next to the points.
    ax_b.annotate(f"n={hi['n']:.0f}/17 — 6 cells die:\n`no downhill origin`",
                  (covs[0], hi['r']), textcoords='offset points', xytext=(18, -26),
                  fontsize=6.6, color=SKETCH_COLOR['FracMinHash'], ha='left')
    ax_b.annotate(f"n={lo['n']:.0f}/17", (covs[0], lo['r']),
                  textcoords='offset points', xytext=(14, -12), fontsize=6.6,
                  color=SKETCH_COLOR['2bRAD anchors'], ha='left')
    ax_b.set_xscale('log')
    ax_b.set_xticks(covs); ax_b.set_xticklabels([f'{c:g}×' for c in covs])
    ax_b.xaxis.set_minor_locator(mticker.NullLocator())
    ax_b.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax_b.set_xlabel('subsampled coverage')
    ax_b.set_title('b  density-matched (s104 vs k16)\nsame estimator, same density',
                   fontsize=8.4, color=INK, pad=5)
    ax_b.legend(loc='lower right', fontsize=7, handletextpad=0.5, labelspacing=0.3)

    # --- c: mechanism — window population at 0.5x ----------------------------
    m = pd.read_csv(f1m, sep='\t')
    m = m[(m['cov'] == 0.5) & (m['arm'].isin(['A_k16', 'E_s104']))]
    rows = {r['arm']: r for _, r in m.iterrows()}
    order = ['A_k16', 'E_s104']
    labs = ['panel k16', 'sketch s104']
    cols = [SKETCH_COLOR['2bRAD anchors'], SKETCH_COLOR['FracMinHash']]
    x = np.arange(2)
    grid(ax_c, axis='y')
    ax_c.bar(x - 0.18, [rows[a]['lm_win_med'] for a in order], width=0.34,
             color=cols, zorder=3)
    ax_c.set_ylabel('landmarks / window (median)')
    ax_c.set_ylim(0, 21)
    for xi, a in zip(x, order):
        ax_c.text(xi - 0.18, rows[a]['lm_win_med'] + 0.5,
                  f"{rows[a]['lm_win_med']:.1f}", ha='center', fontsize=7.5,
                  color=INK)
    ax2 = ax_c.twinx()
    ax2.bar(x + 0.18, [rows[a]['det_frac_med'] for a in order], width=0.34,
            color=cols, alpha=0.42, hatch='//', edgecolor='white', zorder=3)
    ax2.set_ylabel('detected fraction', color=INK2)
    ax2.set_ylim(0, 0.25)
    ax2.tick_params(axis='y', colors=INK2)
    for xi, a in zip(x, order):
        ax2.text(xi + 0.18, rows[a]['det_frac_med'] + 0.008,
                 f"{rows[a]['det_frac_med']:.3f}", ha='center', fontsize=7.5,
                 color=INK2)
    ax_c.set_xticks(x); ax_c.set_xticklabels(labs)
    ax_c.set_title('c  mechanism at 0.5×:\nwindows equally populated',
                   fontsize=8.4, color=INK, pad=5)
    from matplotlib.patches import Patch
    ax_c.legend(handles=[Patch(facecolor=MUTED, label='landmarks / window'),
                         Patch(facecolor=MUTED, alpha=0.42, hatch='//',
                               edgecolor='white', label='detected fraction')],
                loc='upper left', fontsize=6.8, handlelength=1.2)

    ax_a1.text(-0.28, 1.13, 'a', transform=ax_a1.transAxes, fontsize=11,
               fontweight='bold', color=INK, va='top')
    cap = (CAP4 + NL + NL +
           f'(panel a: at 1× the coordinate fit is worth {gain["2bRAD anchors"]:+.2f} r on '
           f'2bRAD anchors but only {gain["FracMinHash"]:+.2f} on a FracMinHash sketch; '
           f'interaction {inter:+.2f}.)')
    fig.text(0.0, -0.02, cap, fontsize=7, color=MUTED, va='top')
    fig.subplots_adjust(left=0.055, right=0.985, top=0.82, bottom=0.14)
    save(fig, 'fig4_attribution')


#: Fig 5 caption: panels e/f keep the GC-sweep text with their new letters.
CAP5 = (
    'Eighteen real genomes spanning GC 25.4\u201372.0%, simulated reads, planted '
    'log\u2082PTR 0.5\u20132.0 (n = 8 cells per point). e, Enzyme-panel landmark density '
    'against GC: k16 rises monotonically (\u22484.5k/Mb at 26% to \u224812.4k/Mb at' + NL +
    '72%) and only the low-GC end is depressed (\u22482\u00d7 below mid-GC); '
    'density-matched FracMinHash tracks the panel within a few percent at every '
    'GC. The k8/k16 ratio thins at the GC extremes (0.63 at 72%) — the extra' + NL +
    'eight enzymes buy high-GC density, not accuracy. f, Pearson r at 1\u00d7: k8 '
    'matches k16 at every GC despite carrying only 63% of its landmarks at 72% '
    'GC (A6 holds across GC), and matched FracMinHash ties the panel.' + NL +
    'Shaded band: the instrument boundary — at GC \u226430% together with 0.5\u00d7 depth '
    'every landmark mode fails (r < 0.7) — so the claimable range is GC \u226530% '
    'at \u22651\u00d7 coverage.'
)


# --- Fig 5: panel design — size vs accuracy/cost, and the panel across GC ---
def fig5():
    """Panel size against accuracy and cost (a–d).

    Depth is an *ordered* variable in a/b/c, so it gets a single-hue sequential
    ramp rather than categorical colours; d compares four *methods*, so it uses
    the fixed categorical slots. e/f: the enzyme panel as an instrument across
    GC — landmark density and 1x accuracy.
    """
    f = ROOT / 'data' / 'panel_sweep.tsv'
    if not f.exists():
        print('  (skipping fig5: panel_sweep.tsv absent)'); return
    d = pd.read_csv(f, sep='\t')
    depths, ks = sorted(d['depth'].unique()), sorted(d['k'].unique())
    ramp = plt.get_cmap('Blues')(np.linspace(0.34, 0.95, len(depths)))

    den = pd.read_csv(ROOT / 'data' / 'f2_gc_sweep' / 'F2_density.tsv', sep='\t')
    acc = pd.read_csv(ROOT / 'data' / 'f2_gc_sweep' / 'F2_accuracy.tsv', sep='\t')
    BLUE, GREEN = SKETCH_COLOR['2bRAD anchors'], SKETCH_COLOR['FracMinHash']
    den = den.copy()
    den['fam'] = den['mode'].map(
        lambda mo: mo if mo.startswith('panel_') else 'fmh_' + mo.rsplit('_', 1)[-1])

    fig, axes = plt.subplots(2, 3, figsize=(10.6, 6.9))
    (ax_r, ax_m, ax_t), (ax_p, ax_gd, ax_ga) = axes
    for ax in axes.ravel():
        grid(ax, axis='both')

    for ax, col, lab, ttl in (
            (ax_r, 'r', 'Pearson r vs growth rate', 'r: tied 4–12, point peak at 8'),
            (ax_m, 'rmse', 'RMSE vs predicted log₂PTR', 'magnitude: best at 4–8 enzymes'),
            (ax_t, 'seconds', 'wall clock per sample (s)', 'cost: linear in k')):
        for c, dp in zip(ramp, depths):
            g = d[d['depth'] == dp].sort_values('k')
            ax.plot(g['k'], g[col], '-o', color=c, ms=4, label=f'{dp:g}×')
        ax.set_xlabel('enzymes in panel'); ax.set_ylabel(lab); ax.set_xticks(ks)
        ax.set_title(ttl, fontsize=8.4, color=INK, pad=4)
    ax_r.legend(title='read depth', title_fontsize=7, ncol=2, loc='lower right',
                handletextpad=0.4, columnspacing=1.0, labelspacing=0.25)
    ax_t.set_ylim(bottom=0)
    # r lives in [0.90, 0.99] here; letting it autoscale would magnify noise that
    # n = 16 media cannot resolve, so the axis is pinned to an honest range.
    ax_r.set_ylim(0.88, 1.0)

    # d: cost against depth, comparing methods rather than panel sizes
    series = [
        (d[d['k'] == 2].sort_values('depth'), 'seconds', ARM_COLOR['A'], 'sk2bGrow, 2 enzymes', '-'),
        (d[d['k'] == 16].sort_values('depth'), 'seconds', ARM_COLOR['B'], 'sk2bGrow, 16 enzymes', '-'),
        (d.drop_duplicates('depth').sort_values('depth'), 'pilea_seconds',
         ARM_COLOR['C_relaxed'], 'Pilea, gates off', '-'),
        (d.drop_duplicates('depth').sort_values('depth'), 'pilea_default_seconds',
         ARM_COLOR['C_default'], 'Pilea, defaults', (0, (4, 2.5))),
    ]
    # Pilea at defaults reports nothing below 10x on this dataset. Plotting its
    # cost there without saying so reads as "5x cheaper" rather than "declined
    # to answer", so those points are hollow — same convention as Fig 2.
    answered = set(grow[(grow['arm'] == 'C_default')].dropna(subset=['log2ptr'])['cov'])
    for g, col, colour, lab, ls in series:
        if col not in g.columns or not np.isfinite(g[col]).any():
            continue
        hollow = 'default' in col
        ax_p.plot(g['depth'], g[col], marker='' if hollow else 'o', ms=4,
                  color=colour, label=lab, ls=ls)
        if hollow:
            for x, y in zip(g['depth'], g[col]):
                filled = x in answered
                ax_p.plot([x], [y], 'o', ms=4.5, color=colour if filled else SURFACE,
                          markeredgecolor=colour, markeredgewidth=1.4, zorder=5)
    ax_p.set_xscale('log'); ax_p.set_xticks(depths)
    ax_p.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f'{v:g}×'))
    ax_p.set_xlabel('read depth'); ax_p.set_ylabel('wall clock per sample (s)')
    ax_p.set_ylim(bottom=0)
    # The amber V fills the frame — ascent, peak and descent all cross every
    # corner a 4-entry legend could use — so the legend and the hollow-note go
    # above the panel, under the title.
    ax_p.legend(loc='lower left', bbox_to_anchor=(0.0, 1.02), fontsize=7,
                handletextpad=0.5, labelspacing=0.25)
    ax_p.text(1.0, 1.02, 'hollow = no estimate returned',
              transform=ax_p.transAxes, fontsize=6.4, color=MUTED,
              va='bottom', ha='right')
    ax_p.set_title('sk2bGrow vs Pilea, same cells\n(single genome; the ordering flips at MAG scale — Fig. 8)',
                   fontsize=8.4, color=INK, pad=56)

    # e: panel density vs GC ---------------------------------------------------
    series = [('panel_k16', 'panel k16', BLUE, '-', 1.0),
              ('panel_k8', 'panel k8', BLUE, '--', 0.55),
              ('panel_k2', 'panel k2', BLUE, ':', 0.35),
              ('fmh_k16', 'FMH matched k16', GREEN, '-', 1.0),
              ('fmh_k8', 'FMH matched k8', GREEN, '--', 0.55)]
    for fam, lab, c, ls, al in series:
        s = den[den['fam'] == fam].sort_values('gc')
        ax_gd.plot(s['gc'], s['landmarks_mb'], ls=ls, color=c, alpha=al, lw=2.0,
                   label=lab, zorder=3)
    k16 = den[den['fam'] == 'panel_k16'].set_index('gc')['landmarks_mb']
    k8 = den[den['fam'] == 'panel_k8'].set_index('gc')['landmarks_mb']
    g72 = max(k16.index)
    ax_gd.annotate(f'k8/k16 = {k8[g72] / k16[g72]:.2f}', xy=(g72, k8[g72]),
                   xytext=(g72 - 24, k8[g72] * 0.55), fontsize=7.5, color=INK2,
                   arrowprops=dict(arrowstyle='-', color=INK2, lw=0.8))
    ax_gd.text(0.03, 0.19, 'low-GC end ≈2× depressed', transform=ax_gd.transAxes,
               fontsize=7, color=MUTED, va='top')
    ax_gd.set_xlabel('genome GC (%)'); ax_gd.set_ylabel('landmarks / Mb')
    ax_gd.set_title('panel density vs GC', loc='left', fontsize=9, color=INK)
    ax_gd.legend(loc='upper left', fontsize=7, handletextpad=0.6, labelspacing=0.3)

    # f: r at 1x, three series -------------------------------------------------
    a = acc[acc['depth'] == 1.0].copy()
    fmh = (a[a['mode'].str.startswith('E_')].groupby('genome')
           .agg(gc=('gc', 'first'), r=('r', 'mean')).reset_index())
    frames = []
    for mode, lab, c, ls in [('A_k16', 'panel k16', BLUE, '-'),
                             ('A_k8', 'panel k8', BLUE, '--')]:
        t = a[a['mode'] == mode][['gc', 'r']].copy()
        t['lab'], t['c'], t['ls'] = lab, c, ls
        frames.append(t)
    fmh['lab'], fmh['c'], fmh['ls'] = 'FMH matched', GREEN, '-'
    frames.append(fmh[['gc', 'r', 'lab', 'c', 'ls']])
    ax_ga.axvspan(24.5, 30.5, color=MUTED, alpha=0.12, zorder=0)
    ax_ga.text(27.5, 0.575, 'dead zone\n(fails at 0.5×)', ha='center',
               fontsize=6.8, color=INK2, va='bottom')
    for t in frames:
        t = t.sort_values('gc')
        ax_ga.plot(t['gc'], t['r'], ls=t['ls'].iloc[0], color=t['c'].iloc[0],
                   lw=1.6, alpha=0.85, zorder=3, label=t['lab'].iloc[0])
        ax_ga.plot(t['gc'], t['r'], 'o', color=t['c'].iloc[0], ms=5, alpha=0.85,
                   markeredgecolor='white', markeredgewidth=0.8, zorder=4)
    ax_ga.set_xlabel('genome GC (%)'); ax_ga.set_ylabel('Pearson r at 1×')
    ax_ga.set_ylim(0.55, 1.02)
    ax_ga.set_title('accuracy at 1× vs GC', loc='left', fontsize=9, color=INK)
    ax_ga.legend(loc='lower right', fontsize=7, handletextpad=0.5, labelspacing=0.3)

    for ax, letter, dx in ((ax_r, 'a', -0.24), (ax_m, 'b', -0.20), (ax_t, 'c', -0.20),
                           (ax_p, 'd', -0.24), (ax_gd, 'e', -0.17), (ax_ga, 'f', -0.17)):
        ax.text(dx, 1.10, letter, transform=ax.transAxes, fontsize=11,
                fontweight='bold', color=INK, va='top')
    fig.text(0.0, -0.02, CAP5, fontsize=7, color=MUTED, va='top')
    fig.subplots_adjust(left=0.09, right=0.985, top=0.945, bottom=0.075,
                        wspace=0.40, hspace=0.65)
    save(fig, 'fig5_panel_design')


#: Fig 6's caption, kept out of the function so the hard line breaks that keep
#: savefig's tight bbox from widening the canvas stay visible.
CAP6 = (
    'The same reads throughout, and 43,707 of the complete genome\'s 43,735 '
    'anchors survive the cut; only the coordinate changes. Fragmenting into 100 '
    'shuffled contigs kills the coordinate fit — every per-enzyme' + NL +
    'V-fit on shuffled windows has r\u00b2 below zero — and the residual estimates are '
    'winner\u2019s-curse artifacts of the origin grid search: permuting window '
    'coordinates and refitting 200\u00d7 gives a null' + NL +
    'distribution of the same order as the estimates on fragmented references, '
    'so the surviving correlation is inflated by the search artifact, not '
    'recovered biology. Slope and RMSE are the primary metrics.' + NL +
    'Every estimate collapses toward zero (panel b, slope {slope:.2f}), and '
    '{qc:.0f}% of those collapsed estimates pass the fusion QC at 5\u201310\u00d7, '
    'because a destroyed coordinate makes every enzyme agree' + NL +
    'there is no gradient. Pilea\'s rank '
    'regression needs no coordinate and is nearly indifferent to the cut, so on '
    'unscaffolded contigs it wins' + NL +
    'outright (r 0.83 against 0.55 at 1×). `sk2bgrow scaffold` is what puts the '
    'coordinate back: it restores the complete-reference curve exactly — the '
    'blue halo under the green — even against' + NL +
    'a different strain. Panel c is why the slope has to be reported beside r: '
    'across 2 to 100 contigs the correlation never leaves 0.86–0.97 while the '
    'slope falls from 0.88 to 0.21. At 50 contigs — a' + NL +
    'draft most people would call good, N50 156 kb — r reads 0.96 and every '
    'estimate is 44% of truth.' + NL +
    'Panel d asks whether a coordinate is needed at all. Fitting the distribution '
    'of window rates instead of their shape in position — log₂ coverage is uniform '
    'across the genome and the width of that' + NL +
    'uniform is log₂(PTR) — recovers most of the loss above 5× (RMSE 0.87 to 0.14) '
    'with no scaffolding, but carries too little information below it. On the '
    'stationary control it reports 0.000 at 1×' + NL +
    'where Pilea reports 1.153, because the per-window standard error is in the '
    'model rather than being read as growth.'
)


# --- Fig 6: does the method survive a fragmented reference? ------------------
def fig6():
    """The MAG case. Panel a is accuracy against coverage per reference
    condition; panel b is estimated against predicted magnitude at 10x, which is
    where the failure shows its shape -- fragmentation does not add noise, it
    collapses every estimate toward zero.

    Colour encodes the reference condition, a different categorical dimension
    from the arms, so it takes its own fixed slot order.
    """
    f = ROOT / 'data' / 'fragmentation.tsv'
    if not f.exists():
        print('  (skipping fig6: fragmentation.tsv absent)'); return
    d = pd.read_csv(f, sep='\t')
    d = d[d['medium'] != 'RUN_OUT']
    covs = sorted(d['cov'].unique())

    fig, axes = plt.subplots(2, 2, figsize=(7.6, 6.4))
    ax_r, ax_m, ax_s, ax_e = axes.ravel()

    grid(ax_r, axis='both')
    for cond in REF_ORDER:
        xs, ys = [], []
        for c in covs:
            s = d[(d['cov'] == c) & (d['cond'] == cond)].dropna(subset=['log2ptr'])
            if len(s) < 3 or s['log2ptr'].nunique() < 2:
                continue
            xs.append(c); ys.append(stats.pearsonr(s['growth_rate'], s['log2ptr'])[0])
        if not xs:
            continue
        # `complete` and `scafRel` coincide almost exactly -- that IS the result,
        # but two curves drawn identically look like one. Draw complete as a wide
        # halo underneath so the agreement is visible rather than hidden.
        wide = cond == 'complete'
        ax_r.plot(xs, ys, 'o-', color=REF_COLOR[cond], label=REF_LABEL[cond],
                  lw=4.0 if wide else 2.0, ms=8 if wide else 5,
                  markeredgecolor='white', markeredgewidth=0.8,
                  zorder=2 if wide else 3)
    ax_r.axhline(0, lw=1.0, color=MUTED, zorder=1)
    ax_r.set_xscale('log')
    ax_r.set_xticks(covs); ax_r.set_xticklabels([f'{c:g}×' for c in covs])
    ax_r.xaxis.set_minor_locator(mticker.NullLocator())
    ax_r.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax_r.set_xlabel('subsampled coverage')
    ax_r.set_ylabel('Pearson r vs measured growth rate\n(r inflated by ori-search '
                    'artifact — read slope / RMSE)', fontsize=8)
    ax_r.set_title('a  ranking', loc='left', fontsize=9, color=INK)
    ax_r.legend(loc='lower right', handletextpad=0.5, labelspacing=0.3)

    grid(ax_m, axis='both')
    top = max(covs)
    lim = (0, 2.05)
    ax_m.plot(lim, lim, ls=':', lw=1.2, color=MUTED, zorder=1)
    ax_m.text(1.98, 1.90, 'y = x', fontsize=7, color=MUTED, ha='right', va='top')
    for cond in REF_ORDER:
        s = d[(d['cov'] == top) & (d['cond'] == cond)].dropna(
            subset=['log2ptr', 'pred_log2ptr'])
        if s.empty:
            continue
        wide = cond == 'complete'
        ax_m.plot(s['pred_log2ptr'], s['log2ptr'], 'o', color=REF_COLOR[cond],
                  ms=9 if wide else 5, markeredgecolor='white',
                  markeredgewidth=0.8, zorder=2 if wide else 3,
                  label=REF_LABEL[cond])
    ax_m.set_xlim(*lim); ax_m.set_ylim(*lim)
    ax_m.set_xlabel('predicted log₂(PTR)')
    ax_m.set_ylabel('estimated log₂(PTR)')
    ax_m.set_title(f'b  magnitude at {top:g}×', loc='left', fontsize=9, color=INK)

    frag = d[(d['cond'] == 'frag') & (d['cov'] >= 5) & np.isfinite(d['log2ptr'])]
    qc = 100 * frag['passed'].mean() if len(frag) else float('nan')
    slope = np.polyfit(*[d[(d['cond'] == 'frag') & (d['cov'] == top)]
                         .dropna(subset=['log2ptr', 'pred_log2ptr'])[c]
                         for c in ('pred_log2ptr', 'log2ptr')], 1)[0]
    fig.text(0.0, -0.02, CAP6.format(qc=qc, slope=slope), fontsize=7,
             color=MUTED, va='top')
    # --- c: where the threshold is ------------------------------------------
    # r and the fitted slope share one 0-1 axis because both are dimensionless
    # and both would be 1 for a perfect estimator. That is the whole point of
    # the panel: they diverge, so r is not a fragmentation diagnostic.
    grid(ax_s, axis='both')
    sw = ROOT / 'data' / 'fragmentation_sweep.tsv'
    if sw.exists():
        w = pd.read_csv(sw, sep='\t')
        w = w[w['medium'] != 'RUN_OUT'].dropna(subset=['log2ptr'])
        ns = sorted(w['n_contigs'].unique())
        rs, sl = [], []
        for n in ns:
            t = w[w['n_contigs'] == n]
            rs.append(stats.pearsonr(t['growth_rate'], t['log2ptr'])[0])
            sl.append(np.polyfit(t['growth_rate'], t['log2ptr'], 1)[0])
        ax_s.axhline(1.0, ls=':', lw=1.0, color=MUTED, zorder=1)
        ax_s.plot(ns, rs, 'o-', color=MUTED, ms=5, markeredgecolor='white',
                  markeredgewidth=0.8, zorder=3)
        ax_s.plot(ns, sl, 'o-', color=REF_COLOR['frag'], ms=5,
                  markeredgecolor='white', markeredgewidth=0.8, zorder=4)
        # The r curve hugs the dotted 1.0 line across the whole frame, so its
        # name cannot sit between curve and line; it goes above the frame at
        # the right end, where nothing else reaches.
        ax_s.annotate('Pearson r — inflated by\nsearch artifact', (ns[-1], rs[-1]),
                      textcoords='offset points', xytext=(-4, 30), ha='right',
                      va='bottom', fontsize=7.5, color=MUTED, annotation_clip=False)
        ax_s.annotate('fitted slope', (ns[2], sl[2]), textcoords='offset points',
                      xytext=(0, -16), ha='center', fontsize=7.5,
                      color=REF_COLOR['frag'])
        ax_s.set_xscale('log')
        ax_s.set_xticks(ns)
        ax_s.set_xticklabels([str(n) for n in ns])
        ax_s.xaxis.set_minor_locator(mticker.NullLocator())
        ax_s.xaxis.set_minor_formatter(mticker.NullFormatter())
        ax_s.set_ylim(0, 1.12)
        ax_s.set_xlabel('contigs the reference is cut into')
        ax_s.set_ylabel('value (1.0 = correct)')
        ax_s.set_title(f'c  correlation cannot see it, at {top:g}×',
                       loc='left', fontsize=9, color=INK, pad=13)

    # --- d: is a coordinate even necessary? ---------------------------------
    # RMSE on a log axis, because the three estimators differ by an order of
    # magnitude and the crossover near 5x is the whole point.
    grid(ax_e, axis='both')
    sf = ROOT / 'data' / 'fragmentation_spread.tsv'
    if sf.exists():
        sp = pd.read_csv(sf, sep='\t')
        e = pd.concat([d[['cond', 'medium', 'cov', 'log2ptr', 'pred_log2ptr']],
                       sp.assign(pred_log2ptr=sp['medium'].map(
                           d.set_index('medium')['pred_log2ptr'].groupby(level=0).first()))
                       [['cond', 'medium', 'cov', 'log2ptr', 'pred_log2ptr']]])
        e = e[e['medium'] != 'RUN_OUT'].dropna(subset=['log2ptr', 'pred_log2ptr'])
        for cond in EST_ORDER:
            xs, ys = [], []
            for c in covs:
                t = e[(e['cov'] == c) & (e['cond'] == cond)]
                if len(t) < 3:
                    continue
                xs.append(c)
                ys.append(float(np.sqrt(((t['log2ptr'] - t['pred_log2ptr']) ** 2).mean())))
            if not xs:
                continue
            ax_e.plot(xs, ys, 'o-', color=EST_COLOR[cond], label=EST_LABEL[cond],
                      ms=5, markeredgecolor='white', markeredgewidth=0.8, zorder=3)
        ax_e.set_xscale('log'); ax_e.set_yscale('log')
        ax_e.set_xticks(covs); ax_e.set_xticklabels([f'{c:g}×' for c in covs])
        for axis in (ax_e.xaxis,):
            axis.set_minor_locator(mticker.NullLocator())
            axis.set_minor_formatter(mticker.NullFormatter())
        ax_e.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f'{v:g}'))
        ax_e.set_xlabel('subsampled coverage')
        ax_e.set_ylabel('RMSE vs predicted log₂(PTR)')
        ax_e.set_title('d  on 100 contigs, which estimator?', loc='left',
                       fontsize=9, color=INK)
        ax_e.legend(loc='lower left', handletextpad=0.5, labelspacing=0.3)

    fig.subplots_adjust(left=0.095, right=0.99, top=0.94, bottom=0.115,
                        wspace=0.28, hspace=0.36)
    save(fig, 'fig6_fragmentation')


#: Fig 7's caption.
CAP7 = (
    'Sun PRJNA689204 stool communities, three samples of \u2248125\u2013135 Gb each. a, '
    'Bland\u2013Altman summary of sk2bGrow-WGS (B) against Pilea at its shipped '
    'defaults (A), on the common denominator of species both methods report' + NL +
    '(n = 58/68/78 per sample): bias is stable at +0.22 to +0.25 log\u2082 (Pilea '
    'reads \u224816\u201319% high) with limits of agreement \u2248 \u22120.43 to +0.89; the low CCC '
    '(0.25\u20130.37) reflects the narrow observed range (span \u22481.4 log\u2082),' + NL +
    'not broken agreement. Per-species points are not committed; each marker is '
    'a sample\u2019s bias with its LoA whisker, annotated with n and the B/A ranges. '
    'b, the real-2bRAD arm (C) against the same WGS arm: exact read-level' + NL +
    'deduplication (95.7% of reads collapse to unique sequences) destroys the '
    'PTR signal — Pearson r falls from 0.73\u20130.80 to 0.30\u20130.60 and CCC collapses '
    '(e.g. S06: r 0.797 to 0.435, CCC 0.728 to 0.192). sk2bGrow\u2019s pipeline must not' + NL +
    'deduplicate before counting. (raw-arm numbers from the c_arm_diag run, '
    'REVIEW \u00a76; deduped from p2_review_agreement.tsv.)'
)


# --- Fig 7: metagenome agreement — BA summary + dedup ablation --------------
def fig7():
    ag = pd.read_csv(ROOT / 'data' / 'sun_three_arm' / 'p2_review_agreement.tsv',
                     sep='\t')
    ba = ag[ag['comparison'] == 'A_default_vs_B'].set_index('sample')
    # raw-arm numbers exist only as prose in REVIEW.md §6 (c_arm_diag, job
    # 4006836); the deduped arms are the committed tsv rows.
    raw = {'S01': (0.449, 0.440), 'S06': (0.797, 0.728), 'S07': (0.734, 0.566)}
    SAMP = [('S01', 'S01'), ('S06', 'S06'), ('S07', 'S07')]
    SC = {'S01': ARM_COLOR['A'], 'S06': ARM_COLOR['B'], 'S07': ARM_COLOR['C_default']}

    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(9.0, 3.5), gridspec_kw={'width_ratios': [1.15, 1]})
    grid(ax_a, axis='y')
    ax_a.axhline(0, lw=1.0, color=MUTED, zorder=1)
    for i, (s, _) in enumerate(SAMP):
        r = ba.loc[s]
        c = SC[s]
        ax_a.plot([i, i], [r['loa_lo'], r['loa_hi']], '-', color=c, lw=1.6,
                  zorder=2)
        ax_a.plot([i, i], [r['loa_lo'], r['loa_hi']], '_', color=c, ms=11,
                  mew=1.6, zorder=2)
        ax_a.plot([i], [r['bias']], 'o', color=c, ms=7, zorder=3,
                  markeredgecolor='white', markeredgewidth=1.0)
        note = f"n={r['n']}\nB [{r['range_x']}]\nA [{r['range_y']}]"
        ax_a.annotate(note, (i, r['bias']), textcoords='offset points',
                      xytext=(-16 if i == 2 else 16, 2),
                      ha='right' if i == 2 else 'left',
                      fontsize=6.6, color=INK2, va='center')
    ax_a.set_xticks(range(3))
    ax_a.set_xticklabels([s for s, _ in SAMP])
    ax_a.set_xlabel('sample')
    ax_a.set_ylabel('\u0394log\u2082(PTR)  (Pilea \u2212 sk2bGrow-WGS)')
    ax_a.set_ylim(-0.75, 1.15)
    ax_a.set_title('method agreement, common denominator', loc='left',
                   fontsize=9, color=INK)
    ax_a.text(0.03, 0.97, 'bias +0.22 to +0.25\nLoA \u2248 \u22120.43 to +0.89',
              transform=ax_a.transAxes, fontsize=6.8, color=MUTED, va='top')

    # b: raw vs deduped, r and CCC
    ded = ag[ag['comparison'] == 'B_vs_C'].set_index('sample')
    grid(ax_b, axis='y')
    x = np.arange(3)
    w = 0.2
    DED = '#bdbdbd'
    for i, (s, _) in enumerate(SAMP):
        rr, cr = raw[s]
        rd, cd = ded.loc[s]['r'], ded.loc[s]['ccc']
        bars = [(rr, SC[s], 0), (rd, DED, 0), (cr, SC[s], 1), (cd, DED, 1)]
        for j, (v, c, k) in enumerate(bars):
            ax_b.bar(i + (j - 1.5) * (w + 0.02), v, width=w, color=c, zorder=3,
                     alpha=0.95 if c == DED else 0.75 + 0.25 * (1 - k))
            ax_b.text(i + (j - 1.5) * (w + 0.02), v + 0.02, f'{v:.2f}',
                      ha='center', fontsize=6.2, color=INK2)
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker='s', ls='', color=SC[s], ms=8, label=s)
               for s, _ in SAMP]
    handles += [Patch(facecolor=DED, label='deduped counts')]
    ax_b.legend(handles=handles, loc='upper right', fontsize=6.8,
                handletextpad=0.4, ncol=2, columnspacing=0.8)
    ax_b.text(0.02, 0.97, 'per sample: left pair = Pearson r,\nright pair = CCC',
              transform=ax_b.transAxes, fontsize=6.8, color=INK2, va='top')
    ax_b.text(0.02, 0.90, 'exact dedup collapses r and CCC',
              transform=ax_b.transAxes, fontsize=6.8, color=MUTED, va='top')
    ax_b.set_xticks(x); ax_b.set_xticklabels([s for s, _ in SAMP])
    ax_b.set_xlabel('sample')
    ax_b.set_ylabel('agreement with WGS arm (B)')
    ax_b.set_ylim(0, 1.05)
    ax_b.set_title('real 2bRAD arm: dedup ablation', loc='left', fontsize=9,
                   color=INK)

    ax_a.text(-0.14, 1.04, 'a', transform=ax_a.transAxes, fontsize=11,
              fontweight='bold', color=INK, va='top')
    ax_b.text(-0.16, 1.04, 'b', transform=ax_b.transAxes, fontsize=11,
              fontweight='bold', color=INK, va='top')
    fig.text(0.0, -0.05, CAP7, fontsize=7, color=MUTED, va='top')
    fig.subplots_adjust(left=0.075, right=0.985, top=0.88, bottom=0.15,
                        wspace=0.34)
    save(fig, 'fig7_metagenome')


#: Fig 8's caption: the MAG-QC caption verbatim, then the cost caption with
#: its panels relettered to c/d.
CAP8 = (
    'RBC metagenome (PRJNA974210): 522 MAGs \u00d7 9 samples. Panels a and '
    'b show the legacy C5 statistics, which predate the current signed '
    'fixed-origin and conservative fragmented-reference policy. a, the '
    'enzyme-consistency QC passes less often as the reference fragments: mean '
    'QC pass rate falls from 22.5% (\u226410 contigs) to 2\u20137% (>25 contigs);' + NL +
    'Spearman \u03c1 = \u22120.41 across 522 MAGs, still \u22120.34 after controlling for '
    'mean coverage (both p < 1e\u22124) — legacy QC tracked fragmentation in '
    'references; coverage is the strongest single predictor and the two' + NL +
    'effects are partially entangled in this MAG population. b, legacy recall under '
    'three denominators: sk2bGrow\u2019s reported fraction is structurally 1.00 '
    '(it has no output gate) — a denominator artifact, not performance; under '
    'the gate-comparable denominator (QC passes / 522) sk2bGrow reports' + NL +
    '3.8\u201313.8% and Pilea (defaults) 5.0\u201311.1%, overlapping ranges; on the '
    'common denominator (the MAGs Pilea reports, which sk2bGrow also '
    'estimates) sk2bGrow\u2019s QC passes 43\u201376% of them.'
    + NL + NL +
    'Cost narrative for scaling past single-genome benchmarks. c, wall clock '
    'per sample relative to Pilea at its defaults: on the 522-MAG RBC dataset '
    'the as-built matcher costs \u224890\u2013240\u00d7 (matching reads against 24.1M' + NL +
    'anchors is 89.6% of count time); dropping to --max-mismatch 1 (two ~16-bp '
    'seeds instead of three ~11-bp seeds) shortens posting lists enough to cut '
    'this to \u22484\u201313\u00d7 measured end-to-end on all nine samples (12\u201328\u00d7 over' + NL +
    'the mm = 2 baseline), at the price of 2.8\u201310.7% of anchors lost (median' + NL +
    '4.4%; 0 genomes lost; unbiased on shared anchors). At GTDB scale (\u226590% of references '
    'absent from a sample) the M4 containment screen adds a 21.1\u00d7 lookup '
    'speedup and drives false positives on absent genomes from 4,048/4,700' + NL +
    'to 0. Under the current auto policy, only 26 of 4,698 '
    'genome\u00d7sample observations pass QC because unsafe sorted-rank '
    'fallback is refused on fragmented MAGs (Table 10). d, GTDB projection '
    '(136,646 species representatives \u2248 629 Gbp, '
    '38.8 B/anchor on disk): enzyme k16 \u2248232 GB, k8 \u2248184 GB, half-density FMH '
    '(scale 200) \u2248124 GB with \u22640.007 r lost at \u22651\u00d7. The often-quoted 752 GB' + NL +
    'is peak resident memory during a build, not disk.'
)


# --- Fig 8: MAG QC — fragmentation filter + recall, and the cost of scaling --
def fig8():
    """RBC metagenome MAG QC (a, b) and the cost narrative for scaling past
    single-genome benchmarks (c, d)."""
    binned = pd.read_csv(ROOT / 'data' / 'c5_review' / 'c5_qc_x_ncontigs_binned.tsv',
                         sep='\t')
    recall = pd.read_csv(ROOT / 'data' / 'c5_review' / 'c5_recall_three_ways.tsv',
                         sep='\t')
    common = pd.read_csv(ROOT / 'data' / 'c5_review' / 'c5_common_denominator.tsv',
                         sep='\t')
    cost = pd.read_csv(ROOT / 'data' / 'c5_review' / 'c5_cost_per_sample.tsv',
                       sep='\t')
    ratio = cost['sk2bgrow_total_min'] / (cost['pilea_default_wall_s'] / 60.0)

    fig, axes = plt.subplots(2, 2, figsize=(7.8, 6.4))
    (ax_a, ax_b), (ax_c, ax_d) = axes
    grid(ax_a, axis='y'); grid(ax_b, axis='y')
    grid(ax_c, axis='y'); grid(ax_d, axis='y')

    # a: QC pass rate vs fragmentation ---------------------------------------
    x = np.arange(len(binned))
    ax_a.bar(x, 100 * binned['qc_rate_mean'], width=0.62,
             color=ARM_COLOR['A'], zorder=3)
    for xi, (_, row) in zip(x, binned.iterrows()):
        ax_a.text(xi, 100 * row['qc_rate_mean'] + 0.6,
                  f"{100 * row['qc_rate_mean']:.0f}%", ha='center', fontsize=7.5,
                  color=INK)
        ax_a.text(xi, -2.8, f"n={row['n_mags']}", ha='center', va='top',
                  fontsize=6.6, color=INK2)
    ax_a.set_xticks(x); ax_a.set_xticklabels(binned['ncontig_bin'], fontsize=7.5)
    ax_a.set_xlabel('contigs in the MAG', labelpad=26)
    ax_a.set_ylabel('QC pass rate (%)')
    ax_a.set_ylim(0, 27)
    ax_a.set_title('legacy QC pass rate vs reference fragmentation', loc='left',
                   fontsize=9, color=INK)
    ax_a.text(0.97, 0.95, 'Spearman \u03c1 = \u22120.41 (raw)\n\u22120.34 controlling coverage\n'
                          '(p < 1e\u22124, n = 522 MAGs)',
              transform=ax_a.transAxes, fontsize=6.6, color=MUTED, va='top',
              ha='right')

    # b: recall under three denominators --------------------------------------
    sk = recall[recall['arm'] == 'sk2bgrow']
    pl = recall[recall['arm'] == 'pilea_default']
    comm = (common['n_common_qc_sk2b'] / common['n_common']).mean()
    groups = [('reported fraction\n(rows output / 522)',
               sk['reported_fraction'].mean(), pl['reported_fraction'].mean()),
              ('gate-comparable\n(QC passes / 522)',
               sk['qc_recall'].mean(), pl['reported_fraction'].mean()),
              ('common denominator\n(QC pass on Pilea\u2019s set)',
               comm, 1.0)]
    x = np.arange(3); w = 0.34
    for k, (lab, sv, pv) in enumerate(groups):
        ax_b.bar(k - w / 2, sv, width=w, color=ARM_COLOR['A'], zorder=3,
                 label='sk2bGrow (legacy)' if k == 0 else None)
        ax_b.bar(k + w / 2, pv, width=w, color=ARM_COLOR['C_default'], zorder=3,
                 label='Pilea (defaults)' if k == 0 else None)
        ax_b.text(k - w / 2, sv + 0.02, f'{sv:.2f}', ha='center', fontsize=7.2,
                  color=ARM_COLOR['A'])
        ax_b.text(k + w / 2, pv + 0.02, f'{pv:.2f}', ha='center', fontsize=7.2,
                  color=ARM_COLOR['C_default'])
    ax_b.set_xticks(x); ax_b.set_xticklabels([g[0] for g in groups], fontsize=7)
    ax_b.set_ylabel('fraction of 522 MAGs')
    ax_b.set_ylim(0, 1.35)
    ax_b.legend(loc='upper right', fontsize=7, handletextpad=0.5)
    ax_b.set_title('legacy recall under three denominators', loc='left', fontsize=9,
                   color=INK)
    ax_b.text(0.30, 0.50, '1.00 = no output gate,\nnot performance', fontsize=6.6,
              color=MUTED)

    # c: cost per sample, staged ----------------------------------------------
    # the mm=1 range is measured end-to-end on all nine samples
    # (data/mm1_e2e/mm1_full_sample_cost.tsv joined to Pilea's per-sample wall)
    mm1 = pd.read_csv(ROOT / 'data' / 'mm1_e2e' / 'mm1_full_sample_cost.tsv',
                      sep='\t')
    pilea_wall = cost[['run', 'pilea_default_wall_s']].rename(
        columns={'run': 'sample'})
    mm1 = mm1.merge(pilea_wall, on='sample')
    mm1_ratio = mm1['wall_s'] / mm1['pilea_default_wall_s']
    stages = [('Pilea\ndefaults', 1.0, 1.0, ARM_COLOR['C_default'], None),
              ('sk2bGrow\nas built\n(522 MAGs)', ratio.min(), ratio.max(),
               ARM_COLOR['A'], None),
              ('+ mm = 1\n(measured)', mm1_ratio.min(), mm1_ratio.max(),
               ARM_COLOR['A'], '//')]
    for i, (lab, lo, hi, c, hatch) in enumerate(stages):
        ax_c.bar(i, hi - lo, bottom=lo, width=0.56, color=c, hatch=hatch,
                 edgecolor='white' if hatch else None, zorder=3)
        top = hi * 1.25
        ax_c.text(i, top, '1×' if hi == lo else f'{lo:.0f}–{hi:.0f}×',
                  ha='center', fontsize=8, color=INK)
    ax_c.set_yscale('log')
    ax_c.set_ylim(0.7, 900)
    ax_c.set_yticks([1, 10, 100])
    ax_c.set_yticklabels(['1', '10', '100'])
    ax_c.set_xticks(range(3)); ax_c.set_xticklabels([s[0] for s in stages],
                                                    fontsize=7.5)
    ax_c.set_ylabel('wall clock per sample (× Pilea defaults)')
    ax_c.set_title('cost per sample, staged', loc='left', fontsize=9, color=INK)
    # The open band between the 1x baseline and the mm=1 bar top (y ~ 20–40)
    # holds the measured-range note; the arrow drops onto the bar's left edge,
    # below the "4–13×" value label it would otherwise cross.
    ax_c.annotate('mm=1 (measured, n=9): 12–28× faster than mm=2 end-to-end;\n'
                  'anchor loss 2.8–10.7% (median 4.4%), 0 genomes lost',
                  xy=(1.72, 8), xytext=(0.02, 0.55),
                  textcoords='axes fraction',
                  fontsize=6.6, color=MUTED, ha='left',
                  arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.8))

    # d: GTDB projection, index disk -------------------------------------------
    disks = [('enzyme k16', 232, ARM_COLOR['A'], None),
             ('enzyme k8', 184, ARM_COLOR['A'], '//'),
             ('FMH half-density\n(scale 200)', 124, SKETCH_COLOR['FracMinHash'],
              None)]
    for i, (lab, gb, c, hatch) in enumerate(disks):
        ax_d.bar(i, gb, width=0.56, color=c, hatch=hatch,
                 edgecolor='white' if hatch else None, zorder=3)
        ax_d.text(i, gb + 6, f'≈{gb} GB', ha='center', fontsize=7.8,
                  color=INK)
    ax_d.set_xticks(range(3)); ax_d.set_xticklabels([d[0] for d in disks],
                                                    fontsize=7.5)
    ax_d.set_ylabel('GTDB index, disk (GB)')
    ax_d.set_ylim(0, 275)
    ax_d.set_title('GTDB projection: index disk', loc='left', fontsize=9,
                   color=INK)
    ax_d.text(0.98, 0.97, 'A1\u2019s 752 GB was peak RSS during the build,\nnot disk '
                          '(disk ≈ 38.8 B/anchor)', transform=ax_d.transAxes,
              fontsize=6.6, color=MUTED, va='top', ha='right')

    for ax, letter, dx in ((ax_a, 'a', -0.14), (ax_b, 'b', -0.13),
                           (ax_c, 'c', -0.12), (ax_d, 'd', -0.14)):
        ax.text(dx, 1.04, letter, transform=ax.transAxes, fontsize=11,
                fontweight='bold', color=INK, va='top')
    fig.text(0.0, -0.02, CAP8, fontsize=7, color=MUTED, va='top')
    fig.subplots_adjust(left=0.10, right=0.985, top=0.93, bottom=0.115,
                        wspace=0.32, hspace=0.42)
    save(fig, 'fig8_mag_qc_cost')


if __name__ == '__main__':
    print('regenerating figures ->', OUT)
    fig1(); fig2(); fig3(); fig4(); fig5(); fig6(); fig7(); fig8()
    print('done')
