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


if __name__ == '__main__':
    print('generating tables ->', OUT)
    table1(); table2(); table3()
    print('done')
