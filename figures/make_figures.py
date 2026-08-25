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
import matplotlib.ticker as mticker
import matplotlib.patheffects as pe

sys.path.insert(0, str(Path(__file__).parent))
from style import (ARM_COLOR, ARM_LABEL, ARM_ORDER, ATTRIBUTION, INK, INK2, MUTED,
                   SKETCH_COLOR, GRID as GRIDC, SURFACE, apply, grid)

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


# --- Fig 2: accuracy against measured growth rate, across coverage -----------
def fig2():
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
    ax.xaxis.set_minor_locator(mticker.NullLocator())
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax.xaxis.set_minor_locator(mticker.NullLocator())
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax.set_xlabel('subsampled coverage')
    ax.set_ylabel('Pearson r vs measured growth rate')
    # arm B collapses to r = 0.16 at 0.5x; the axis has to reach it, or the
    # series would leave the frame with no indication that it did.
    ax.set_ylim(0.1, 1.03); ax.set_xlim(0.42, 13)
    ax.set_title('Accuracy against an independent ground truth', loc='left', color=INK)
    ax.legend(loc='lower right')
    fig.text(0.0, -0.05,
             'Zheng et al. 2020 E. coli K-12, 16 growth media (λ = 0.40–1.72 h⁻¹). '
             'Pilea at its shipped defaults\nreturns no estimate below 10×; at 0.5× with gates off it '
             'returns PTR = 1.0 for every sample, so\nneither has a defined correlation there.',
             fontsize=7, color=MUTED, va='top')
    save(fig, 'fig2_accuracy_vs_coverage')


# --- Fig 3: is the magnitude right, not just the ranking? --------------------
def fig3():
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
    save(fig, 'fig3_magnitude')


# --- Fig 4: the negative control --------------------------------------------
def fig4():
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
    save(fig, 'fig4_negative_control')


# --- Fig 5: what is responsible — sketch or estimator? ----------------------
def fig5():
    """The full 2x2. One panel per estimator, one line per sketch, so the
    interaction is the difference between the two panels rather than something
    the reader has to compute.

    Colour encodes the *sketch* and is identical across panels: the same entity
    keeps the same hue, and the panel carries the estimator.
    """
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.2), sharey=True)
    covs = sorted(grow['cov'].unique())
    ESTIMATORS = [('coordinate V-fit', 'coordinate V-fit (ours)'),
                  ('rank regression', 'rank regression (iRep / Pilea)')]

    def r_at(arm, c):
        s_ = grow[(grow['cov'] == c) & (grow['arm'] == arm)].dropna(subset=['log2ptr'])
        if len(s_) < 3 or s_['log2ptr'].nunique() < 2:
            return np.nan
        return stats.pearsonr(s_['growth_rate'], s_['log2ptr'])[0]

    for ax, (est, title) in zip(axes, ESTIMATORS):
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
        # apart -- otherwise they collide exactly where the lines do.
        shared = [c for c in covs if all(c in d for d in drawn.values())]
        if shared and len(drawn) == 2:
            d1, d2 = drawn.values()
            cx = max(shared, key=lambda c: abs(d1[c] - d2[c]))
            first = cx == covs[0]
            nxt = [c for c in covs if c > cx]
            for sketch, d in drawn.items():
                # Sit the label on the side the curve is leaving, so a rising
                # series never has its own line drawn through its name.
                rising = bool(nxt) and d.get(nxt[0], d[cx]) > d[cx]
                ax.annotate(sketch, (cx, d[cx]), textcoords='offset points',
                            xytext=(3 if first else 0, -15 if rising else 9),
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
    axes[0].set_ylabel('Pearson r vs measured growth rate')
    axes[0].set_ylim(0.1, 1.04)

    at1 = {k: r_at(v, 1.0) for k, v in ATTRIBUTION.items()}
    gain = {sk: at1[('coordinate V-fit', sk)] - at1[('rank regression', sk)]
            for sk in SKETCH_COLOR}
    inter = gain['2bRAD anchors'] - gain['FracMinHash']
    cap = (
        f"The two factors are not additive. At 1x the coordinate fit is worth "
        f"{gain['2bRAD anchors']:+.2f} r on 2bRAD anchors but only "
        f"{gain['FracMinHash']:+.2f} on a FracMinHash sketch" + NL +
        f'(interaction {inter:+.2f}). The sketch effect changes sign with the estimator: anchors are '
        'ahead under the coordinate fit and' + NL +
        'behind under rank regression. Neither component carries the result on its own.'
    ).replace('1x', '1\u00d7')
    fig.text(0.0, -0.02, cap, fontsize=7, color=MUTED, va='top')
    fig.subplots_adjust(left=0.095, right=0.985, top=0.88, bottom=0.17, wspace=0.09)
    save(fig, 'fig5_attribution')


# --- Fig 6: multi-strain simulation — accuracy, recall and cost -------------
def fig6():
    """Three panels because three different things matter and they trade off:
    how often a method answers, how right it is, and what it costs."""
    f = ROOT / 'data' / 'sim_results.tsv'
    if not f.exists():
        print('  (skipping fig6: sim_results.tsv absent)'); return
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
    axes[0].legend(loc='lower right')
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
    save(fig, 'fig6_simulation')


# --- Fig 7: how many enzymes does the panel actually need? ------------------
def fig7():
    """Panel size against accuracy and cost.

    Depth is an *ordered* variable in a/b/c, so it gets a single-hue sequential
    ramp rather than categorical colours; d compares four *methods*, so it uses
    the fixed categorical slots.
    """
    f = ROOT / 'data' / 'panel_sweep.tsv'
    if not f.exists():
        print('  (skipping fig7: panel_sweep.tsv absent)'); return
    d = pd.read_csv(f, sep='\t')
    depths, ks = sorted(d['depth'].unique()), sorted(d['k'].unique())
    ramp = plt.get_cmap('Blues')(np.linspace(0.34, 0.95, len(depths)))

    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.2))
    (ax_r, ax_e), (ax_t, ax_p) = axes
    for ax in axes.ravel():
        grid(ax, axis='both')

    for ax, col, lab, ttl in (
            (ax_r, 'r', 'Pearson r vs growth rate', 'ranking: peaks at 8 enzymes'),
            (ax_e, 'rmse', 'RMSE vs measured log₂PTR', 'magnitude: best at 4–8 enzymes'),
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
            ax_p.text(0.98, 0.05, 'hollow = no estimate returned',
                      transform=ax_p.transAxes, fontsize=6.4, color=MUTED, ha='right')
    ax_p.set_xscale('log'); ax_p.set_xticks(depths)
    ax_p.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f'{v:g}×'))
    ax_p.set_xlabel('read depth'); ax_p.set_ylabel('wall clock per sample (s)')
    ax_p.set_ylim(bottom=0)
    ax_p.legend(loc='upper left', handletextpad=0.5, labelspacing=0.25)
    ax_p.set_title('sk2bGrow vs Pilea, same cells', fontsize=8.4, color=INK, pad=4)

    for ax, letter in zip(axes.ravel(), 'abcd'):
        ax.text(-0.22, 1.11, letter, transform=ax.transAxes, fontsize=11,
                fontweight='bold', color=INK, va='top')
    fig.subplots_adjust(left=0.10, right=0.985, top=0.93, bottom=0.095,
                        wspace=0.33, hspace=0.46)
    save(fig, 'fig7_panel_size')


if __name__ == '__main__':
    print('regenerating figures ->', OUT)
    fig1(); fig2(); fig3(); fig4(); fig5(); fig6(); fig7()
    print('done')
