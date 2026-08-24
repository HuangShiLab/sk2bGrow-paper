#!/usr/bin/env python3
"""Generate the manuscript tables as markdown + TSV from the committed data.

    python3 figures/make_tables.py

Same contract as make_figures.py: reads only data/*.tsv, writes tables/.
"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA, OUT = ROOT / 'data', ROOT / 'tables'
OUT.mkdir(exist_ok=True)


def write(name, df, title, note=''):
    df.to_csv(OUT / f'{name}.tsv', sep='\t', index=False, float_format='%.4f')
    with open(OUT / f'{name}.md', 'w') as fh:
        fh.write(f'**{title}**\n\n{df.to_markdown(index=False, floatfmt=".3f")}\n')
        if note:
            fh.write(f'\n{note}\n')
    print(f'  {name}: {len(df)} rows')


# --- Table 1: enzyme panel, measured vs the design report -------------------
def table1():
    f = DATA / 'density_3genomes.tsv'
    if not f.exists():
        print('  (skipping table1: density_3genomes.tsv absent)'); return
    d = pd.read_csv(f, sep='\t')
    write('table1_enzyme_panel', d,
          'Table 1. Anchor density of the 16-enzyme panel across three genomes',
          'Densities in tags per Mb. Measured with the panel transcribed from '
          '2bRADExtraction.pl. 46 of 48 cells reproduce the design report within 3%; '
          'HaeIV differs by exactly 2.00x (locus vs window counting) and Hin4I is '
          'unreconciled.')


# --- Table 2: accuracy on the E. coli isolate dataset -----------------------
def table2():
    res = pd.read_csv(DATA / 'results_raw.tsv', sep='\t')
    grow = res[res['medium'] != 'RUN_OUT']
    ARM = {'A': 'sk2bGrow', 'B': 'sk2bGrow (Pilea-parity estimator)',
           'C_default': 'Pilea (defaults)', 'C_relaxed': 'Pilea (gates off)'}
    rows = []
    for cov in sorted(grow['cov'].unique()):
        for arm, lab in ARM.items():
            s = grow[(grow['cov'] == cov) & (grow['arm'] == arm)].dropna(subset=['log2ptr'])
            if len(s) < 3 or s['log2ptr'].nunique() < 2:
                rows.append(dict(coverage=cov, method=lab, n=len(s), pearson_r=np.nan,
                                 rmse_vs_predicted=np.nan, slope=np.nan))
                continue
            ok = s.dropna(subset=['pred_log2ptr'])
            rows.append(dict(
                coverage=cov, method=lab, n=len(s),
                pearson_r=stats.pearsonr(s['growth_rate'], s['log2ptr'])[0],
                rmse_vs_predicted=float(np.sqrt(((ok['log2ptr'] - ok['pred_log2ptr']) ** 2).mean())),
                slope=float(np.polyfit(s['growth_rate'], s['log2ptr'], 1)[0])))
    write('table2_ecoli_accuracy', pd.DataFrame(rows),
          'Table 2. Accuracy on the Zheng et al. E. coli dataset, by coverage',
          'Pearson r is against independently measured growth rate. RMSE is against '
          'the lambda*C-derived prediction, which is NOT independent (derived from the '
          'same reads by marker-frequency analysis). Blank rows: the method returned '
          'no estimate, or a constant.')


# --- Table 3: multi-strain simulation, accuracy + cost ----------------------
def table3():
    f = DATA / 'sim_results.tsv'
    if not f.exists():
        print('  (skipping table3: sim_results.tsv absent)'); return
    d = pd.read_csv(f, sep='\t')
    g = (d.groupby(['arm', 'n_strains', 'coverage'])
           .agg(recall=('recall', 'mean'), rmse=('rmse', 'mean'), bias=('bias', 'mean'),
                seconds=('seconds', 'mean'), peak_rss_mb=('peak_rss_mb', 'mean'))
           .reset_index())
    write('table3_simulation', g,
          'Table 3. Multi-strain simulation: accuracy and computational cost',
          '**recall** is the fraction of truly-present strains for which the method '
          'returned any estimate, and must be read alongside RMSE — a method that '
          'reports only the easy cases earns a flattering RMSE. Pilea at its shipped '
          'defaults returns nothing below 8x.')

    agg = (d.groupby('arm')
             .agg(recall=('recall', 'mean'), rmse=('rmse', 'mean'), bias=('bias', 'mean'),
                  spurious=('spurious', 'mean'), seconds=('seconds', 'mean'),
                  peak_rss_mb=('peak_rss_mb', 'mean')).reset_index())
    write('table4_simulation_aggregate', agg,
          'Table 4. Multi-strain simulation, aggregate over the whole grid',
          '**spurious** counts genomes reported that were not in the sample (false '
          'positives); both methods scored zero.')


# --- Table 5: which Pilea gate actually suppresses the estimate -------------
def table5():
    """Pilea's gates-off run still reports the three gated statistics per genome,
    so re-applying each default threshold to that output is an exact ablation and
    needs no extra Pilea runs."""
    f = DATA / 'pilea_gate_statistics.tsv'
    if not f.exists():
        print('  (skipping table5: pilea_gate_statistics.tsv absent)'); return
    d = pd.read_csv(f, sep='\t')
    GATES = [('min-cove (-x 5)', 'coverage', 5.0),
             ('min-frac (-z 0.75)', 'fraction', 0.75),
             ('min-cont (-c 0.25)', 'containment', 0.25)]
    for label, col, thr in GATES:
        d[label] = d[col] >= thr
    d['all gates'] = np.logical_and.reduce([d[l] for l, _, _ in GATES])

    rows = []
    for (ds, depth), g in d.groupby(['dataset', 'depth']):
        rec = {'dataset': ds, 'depth': depth, 'n': len(g),
               'median k-mer coverage': g['coverage'].median()}
        rec.update({l: 100 * g[l].mean() for l, _, _ in GATES})
        rec['reported by defaults'] = 100 * g['all gates'].mean()
        rows.append(rec)
    t = pd.DataFrame(rows).sort_values(['dataset', 'depth'])
    fail = d[~d['all gates']]
    sole = {l: int(((~fail[l]) & np.logical_and.reduce(
                [fail[o] for o, _, _ in GATES if o != l])).sum()) for l, _, _ in GATES}
    write('table5_pilea_gate_ablation', t,
          'Table 5. Which of Pilea\'s quality gates suppresses the estimate',
          'Percentages are the share of genome-estimates passing that gate alone; '
          '"reported by defaults" is the share passing all three. Of the '
          f'{len(fail)} estimates the defaults discard, the sole cause is '
          + ', '.join(f'{l} for {n}' for l, n in sole.items()) + '. '
          'A 150 bp read yields 120 31-mers, so the k-mer-coverage threshold of 5 '
          'corresponds to about 6.6x read coverage.')


# --- Table 6: single-enzyme accuracy, and the ranking the sweep uses --------
def table6():
    """Every sk2bGrow run writes an independent V-shape fit per enzyme before
    fusion, so the single-enzyme benchmark is a re-read of runs already done.

    Caveat stated in the note: these fits share an origin estimated from all 16
    enzymes pooled, so the column answers "how informative is this enzyme given a
    good origin", not "how would a panel of only this enzyme behave". Table 7
    tests the latter directly."""
    f = DATA / 'per_enzyme_zheng.tsv'
    if not f.exists():
        print('  (skipping table6: per_enzyme_zheng.tsv absent)'); return
    d = pd.read_csv(f, sep='\t')
    grow, ctl = d[d.medium != 'RUN_OUT'], d[d.medium == 'RUN_OUT']
    depths = sorted(d.depth.unique())

    rows = []
    for e in sorted(d.enzyme.unique()):
        se = d[d.enzyme == e]
        rec = {'enzyme': e, 'anchors': int(se.na.median()),
               'windows': int(se.nw.median()), 'fit_rate': se.ok.mean()}
        for dp in depths:
            g = grow[(grow.enzyme == e) & (grow.depth == dp) & np.isfinite(grow.log2ptr)]
            rec[f'r@{dp:g}x'] = stats.pearsonr(g.growth_rate, g.log2ptr)[0] if len(g) > 2 else np.nan
            rec[f'rmse@{dp:g}x'] = float(np.sqrt(np.mean((g.log2ptr - g.pred) ** 2))) if len(g) > 2 else np.nan
        c = ctl[ctl.enzyme == e]
        rec['ctl_bias'] = float(np.nanmean(np.abs(c.log2ptr))) if len(c) else np.nan
        low = [dp for dp in depths if dp <= 2]
        rec['r_low'] = float(np.nanmean([rec[f'r@{dp:g}x'] for dp in low]))
        rec['rmse_low'] = float(np.nanmean([rec[f'rmse@{dp:g}x'] for dp in low]))
        rows.append(rec)
    t = pd.DataFrame(rows)
    # low-depth accuracy is what the panel is for; the control bias and the
    # magnitude error are penalties, not tie-breakers on correlation alone
    t['score'] = t.r_low - 0.25 * t.rmse_low - 0.25 * t.ctl_bias
    t = t.sort_values('score', ascending=False).reset_index(drop=True)
    t.insert(0, 'rank', t.index + 1)
    cols = ['rank', 'enzyme', 'anchors', 'windows', 'fit_rate'] + \
           [f'r@{dp:g}x' for dp in depths] + ['r_low', 'rmse_low', 'ctl_bias', 'score']
    write('table6_enzyme_ranking', t[cols],
          'Table 6. Single-enzyme PTR accuracy on the E. coli panel',
          'Each enzyme fitted alone, before fusion; r is against measured growth rate '
          'across 16 media. r_low and rmse_low average the 0.5/1/2x depths; ctl_bias is '
          'mean |log2PTR| on the run-out control, which should be 0. The ranking tracks '
          'anchor yield almost exactly: the top six are the six enzymes with more than '
          '2,500 anchors, the bottom four have fewer than 450 and are anti-correlated at '
          '0.5x. These fits share an origin estimated from all 16 enzymes pooled, so the '
          'table measures how informative an enzyme is given a good origin, not how a '
          'panel of that enzyme alone would behave -- Table 7 tests that directly.')


if __name__ == '__main__':
    print('generating tables ->', OUT)
    table1(); table2(); table3(); table5(); table6()
    print('done')
