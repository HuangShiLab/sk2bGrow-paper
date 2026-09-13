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
    # tabulate applies a scalar floatfmt to integer columns too, which turns a
    # count like 17055 into 17055.000. Give it one format per column instead.
    def fmt(c):
        if pd.api.types.is_integer_dtype(df[c]):
            return '.0f'
        # coverage columns hold 0.5 -- one decimal reads as a depth, three
        # reads as a measurement it is not
        return '.1f' if c in ('coverage', 'cov', 'depth_x') else '.3f'

    fmts = [fmt(c) for c in df.columns]
    with open(OUT / f'{name}.md', 'w') as fh:
        fh.write(f'**{title}**\n\n{df.to_markdown(index=False, floatfmt=fmts)}\n')
        if note:
            fh.write(f'\n{note}\n')
    print(f'  {name}: {len(df)} rows')


def write_blocks(name, blocks, title, note=''):
    """Multi-block variant of write(): one table per block in the .md (each
    under a bold subtitle, skipped when the subtitle is empty), the .tsv
    stacks the same frames separated by '# block:' comment lines."""
    md, tsv = [f'**{title}**\n'], []
    for subtitle, df in blocks:
        def fmt(c):
            if pd.api.types.is_integer_dtype(df[c]):
                return '.0f'
            return '.1f' if c in ('coverage', 'cov', 'depth_x') else '.3f'

        fmts = [fmt(c) for c in df.columns]
        if subtitle:
            md.append(f'**{subtitle}**\n')
        md.append(df.to_markdown(index=False, floatfmt=fmts))
        md.append('')
        tsv.append(f'# block: {subtitle or "main"}')
        # no float_format: pandas' default shortest-round-trip rendering keeps
        # regenerated files byte-stable (the T1 audit block predates %.4f)
        tsv.append(df.to_csv(sep='\t', index=False).rstrip('\n'))
        tsv.append('')
    with open(OUT / f'{name}.md', 'w') as fh:
        fh.write('\n'.join(md).rstrip() + '\n')
        if note:
            fh.write(f'\n{note}\n')
    with open(OUT / f'{name}.tsv', 'w') as fh:
        fh.write('\n'.join(tsv).rstrip() + '\n')
    print(f'  {name}: {sum(len(df) for _, df in blocks)} rows '
          f'in {len(blocks)} blocks')


# --- Table 1: enzyme panel, measured vs the design report -------------------
def table1():
    f = DATA / 'density_3genomes.tsv'
    if not f.exists():
        print('  (skipping table1: density_3genomes.tsv absent)'); return
    audit = pd.read_csv(f, sep='\t')
    blocks = [('', audit,)]
    note = ('Densities in tags per Mb. Measured with the panel transcribed from '
            '2bRADExtraction.pl. 46 of 48 cells reproduce the design report within 3%; '
            'HaeIV differs by exactly 2.00x (locus vs window counting) and Hin4I is '
            'unreconciled. ')

    # Second block: the F2 GC sweep. One row per GC level (median across the
    # genomes at that level), panel subsets k2..k16 plus the density-matched
    # FracMinHash arms.
    f2 = DATA / 'f2_gc_sweep' / 'F2_density.tsv'
    if f2.exists():
        d = pd.read_csv(f2, sep='\t')
        # canonical GC levels as reported in the F2 REVIEW (near-identical
        # genomes merged into one row)
        CANON = [25.4, 26.1, 30.6, 33.0, 33.5, 43.5, 50.8, 57.0, 61.5,
                 65.3, 66.0, 72.0]
        d['gc_level'] = d['gc'].map(lambda g: min(CANON, key=lambda c: abs(c - g)))
        rows = []
        for gc, g in d.groupby('gc_level'):
            rec = {'gc_pct': gc, 'n_genomes': g['genome'].nunique()}
            for mode, col in (('panel_k2', 'panel_k2'), ('panel_k4', 'panel_k4'),
                              ('panel_k8', 'panel_k8'), ('panel_k16', 'panel_k16')):
                rec[col] = g.loc[g['mode'] == mode, 'landmarks_mb'].median()
            rec['k8_k16'] = rec['panel_k8'] / rec['panel_k16']
            rec['fmh_k8'] = g.loc[g['mode'].str.endswith('match_k8'), 'landmarks_mb'].median()
            rec['fmh_k16'] = g.loc[g['mode'].str.endswith('match_k16'), 'landmarks_mb'].median()
            rows.append(rec)
        gc_block = pd.DataFrame(rows).sort_values('gc_pct').reset_index(drop=True)
        blocks.append(('GC sweep: panel density across 18 genomes (GC 25.4–72.0%)',
                       gc_block))
        note += ('GC-sweep block: panel k16 density rises with GC from ~4.5k/Mb at '
                 '26% to ~12.4k/Mb at 72% (low-GC depression only, no high-GC '
                 'collapse); the k8/k16 ratio thins at the GC extremes (0.72 at '
                 '30.6%, 0.64 at 72%). Density-matched FracMinHash arms track the '
                 'panel within 4% at every GC level. Genomes are GTDB R232 '
                 'stand-ins (R226 metadata unavailable locally); reads are '
                 'self-simulated with a planted replication gradient (F2 REVIEW '
                 'protocol); each cell is the median of the genomes sharing a GC '
                 'level.')
    write_blocks('table1_enzyme_panel', blocks,
                 'Table 1. Anchor density of the 16-enzyme panel across three genomes',
                 note)


# --- Table 2: accuracy on the E. coli isolate dataset -----------------------
def table2():
    res = pd.read_csv(DATA / 'results_raw.tsv', sep='\t')
    grow = res[res['medium'] != 'RUN_OUT']
    # Ordered so the rows read as the 2x2 of Methods S4: both estimators on
    # both sketches, then Pilea at its shipped gates.
    ARM = {'A': 'sk2bGrow: anchors + V-fit',
           'E': 'FracMinHash + V-fit',
           'B': 'anchors + rank regression',
           'C_relaxed': 'Pilea gates off: FracMinHash + rank',
           'C_default': 'Pilea (defaults)'}
    rows = []
    for cov in sorted(grow['cov'].unique()):
        for arm, lab in ARM.items():
            s = grow[(grow['cov'] == cov) & (grow['arm'] == arm)].dropna(subset=['log2ptr'])
            if len(s) < 3 or s['log2ptr'].nunique() < 2:
                rows.append(dict(coverage=cov, method=lab, n=len(s), pearson_r=np.nan,
                                 rmse_vs_predicted=np.nan, bias=np.nan, slope=np.nan))
                continue
            ok = s.dropna(subset=['pred_log2ptr'])
            rows.append(dict(
                coverage=cov, method=lab, n=len(s),
                pearson_r=stats.pearsonr(s['growth_rate'], s['log2ptr'])[0],
                rmse_vs_predicted=float(np.sqrt(((ok['log2ptr'] - ok['pred_log2ptr']) ** 2).mean())),
                bias=float((ok['log2ptr'] - ok['pred_log2ptr']).mean()),
                slope=float(np.polyfit(s['growth_rate'], s['log2ptr'], 1)[0])))
    ci = pd.read_csv(DATA / 'ci_bootstrap' / 'table2_ci.tsv', sep='\t')
    ci = ci[ci['arm'].isin(ARM)].rename(columns={'cov': 'coverage'})
    ci['method'] = ci['arm'].map(ARM)
    ci = ci[['coverage', 'method', 'r_lo', 'r_hi', 'slope_lo', 'slope_hi']]
    df = pd.DataFrame(rows).merge(ci, on=['coverage', 'method'], how='left')
    write('table2_ecoli_accuracy', df,
          'Table 2. Accuracy on the Zheng et al. E. coli dataset, by coverage',
          'Pearson r is against independently measured growth rate. RMSE is against '
          'the lambda*C-derived prediction, which is NOT independent (derived from the '
          'same reads by marker-frequency analysis). Blank rows: the method returned '
          'no estimate, or a constant. The first four methods are the 2x2 of sketch '
          '(2bRAD anchors vs FracMinHash) by estimator (coordinate V-fit vs sorted-rank '
          'regression), all run on the same subsampled reads. pearson_r and slope '
          'carry 95% CIs from a media bootstrap (10,000 resamples of the 16 media; '
          'data/ci_bootstrap/table2_ci.tsv).')


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
          'corresponds to about 6.25x read coverage (5 x 150/120).')


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
          'anchor yield with exceptions: the top six by accuracy all carry >1,900 anchors '
          '(BslFI, 2,501 anchors, ranks ninth); the bottom four by accuracy (FalI/AloI/'
          'BplI/PsrI, 687/441/379/419 anchors) are near-zero or negative at 0.5x; PpiI '
          '(336 anchors) is the exception at r 0.61 at 0.5x. These fits share an origin estimated from all 16 enzymes pooled, so the '
          'table measures how informative an enzyme is given a good origin, not how a '
          'panel of that enzyme alone would behave -- Table 7 tests that directly.')


# --- Table 7: panel size vs accuracy and cost -------------------------------
def table7():
    f = DATA / 'panel_sweep.tsv'
    if not f.exists():
        print('  (skipping table7: panel_sweep.tsv absent)'); return
    d = pd.read_csv(f, sep='\t')
    low = d[d['depth'] <= 2]
    g = low.groupby('k').agg(anchors=('anchors', 'first'), r=('r', 'mean'),
                             RMSE=('rmse', 'mean'), slope=('slope', 'mean'),
                             ctl_bias=('ctl_bias', 'mean')).reset_index()
    allk = d.groupby('k').agg(index_s=('index_seconds', 'first'),
                              profile_s=('seconds', 'mean'),
                              peak_RSS_MB=('rss_mb', 'mean')).reset_index()
    t = g.merge(allk, on='k')
    t['speedup_vs_16'] = t.loc[t['k'] == 16, 'profile_s'].iloc[0] / t['profile_s']
    for c in ('k', 'anchors'):
        t[c] = t[c].astype(int)
    t['peak_RSS_MB'] = t['peak_RSS_MB'].round().astype(int)
    t = t.rename(columns={'k': 'enzymes', 'r': 'r (<=2x)', 'RMSE': 'RMSE (<=2x)',
                          'slope': 'slope (<=2x)', 'ctl_bias': 'run-out bias (<=2x)'})
    ci7 = pd.read_csv(DATA / 'ci_bootstrap' / 'panel_sweep_ci.tsv', sep='\t')
    ci7 = ci7[pd.to_numeric(ci7['k'], errors='coerce').notna()].copy()
    ci7['enzymes'] = ci7['k'].astype(int)
    t = t.merge(ci7[['enzymes', 'r_lo', 'r_hi']], on='enzymes', how='left')
    t = t.rename(columns={'r_lo': 'r_lo (<=2x)', 'r_hi': 'r_hi (<=2x)'})

    pil = []
    for lab, sc, rc in (('Pilea, gates off', 'pilea_seconds', 'pilea_rss_mb'),
                        ('Pilea, defaults', 'pilea_default_seconds', 'pilea_default_rss_mb')):
        if sc in d.columns and np.isfinite(d[sc]).any():
            u = d.drop_duplicates('depth')
            pil.append(f'{lab} {np.nanmean(u[sc]):.1f} s / {np.nanmean(u[rc]):.0f} MB')
    write('table7_panel_size', t,
          'Table 7. Enzyme-panel size against accuracy and cost (E. coli, 16 media)',
          'Subsets are the top k of Table 6. Accuracy columns average the 0.5/1/2x depths, '
          'the regime the panel exists for; cost columns average all five. Index cost is '
          'one-off per reference. Accuracy is a plateau: 4-12 enzymes are statistically '
          'indistinguishable (bootstrap 95% CIs overlap; k = 8 is the point-estimate '
          'peak), while '
          'cost is linear in k -- the four sparsest enzymes (< 450 anchors) buy nothing and '
          'double the bias on the replication run-out control. We did not establish the '
          'mechanism: two candidate explanations were tested and both refuted (sparse '
          'enzymes are not more biased than dense ones, and forcing fixed-effect weights '
          'does not recover the loss). '
          + ('Same cells, ' + '; '.join(pil) + '.' if pil else ''))


# --- Table 8: reference fragmentation ---------------------------------------
def table8():
    f = DATA / 'fragmentation.tsv'
    if not f.exists():
        print('  (skipping table8: fragmentation.tsv absent)'); return
    d = pd.read_csv(f, sep='\t')
    sf = DATA / 'fragmentation_spread.tsv'
    if sf.exists():
        sp = pd.read_csv(sf, sep='\t')
        pred = d.drop_duplicates('medium').set_index('medium')['pred_log2ptr']
        gr = d.drop_duplicates('medium').set_index('medium')['growth_rate']
        sp['pred_log2ptr'] = sp['medium'].map(pred)
        sp['growth_rate'] = sp['medium'].map(gr)
        d = pd.concat([d, sp], ignore_index=True)

    LAB = {'complete': 'complete chromosome',
           'frag': '100 contigs',
           'scafSelf': '100 contigs, scaffolded vs itself',
           'scafRel': '100 contigs, scaffolded vs a relative',
           'pileaFrag': '100 contigs, Pilea',
           'spread_frag': '100 contigs, order-free spread MLE'}
    grow = d[d['medium'] != 'RUN_OUT']
    ctl = d[d['medium'] == 'RUN_OUT'].set_index(['cond', 'cov'])['log2ptr']
    rows = []
    for cov in sorted(grow['cov'].unique()):
        for cond, lab in LAB.items():
            s = grow[(grow['cov'] == cov) & (grow['cond'] == cond)].dropna(
                subset=['log2ptr'])
            if s.empty:
                continue
            ok = s.dropna(subset=['pred_log2ptr'])
            e = ok['log2ptr'] - ok['pred_log2ptr']
            degenerate = len(s) < 3 or s['log2ptr'].nunique() < 2
            rows.append(dict(
                coverage=cov, reference=lab, n=len(s),
                pearson_r=np.nan if degenerate else
                stats.pearsonr(s['growth_rate'], s['log2ptr'])[0],
                rmse=float(np.sqrt((e ** 2).mean())) if len(ok) else np.nan,
                bias=float(e.mean()) if len(ok) else np.nan,
                slope=np.nan if degenerate else
                float(np.polyfit(s['growth_rate'], s['log2ptr'], 1)[0]),
                # Only sk2bGrow has a fusion QC. analyze.py defaults `passed`
                # to True where the column is absent, which would print a
                # spurious 100% for Pilea in the one table that is about the QC.
                qc_pass=100 * s['passed'].mean()
                if cond.startswith(('complete', 'frag', 'scaf'))
                and 'passed' in s and s['passed'].notna().any() else np.nan,
                stationary_control=ctl.get((cond, cov), np.nan)))
    write('table8_fragmentation', pd.DataFrame(rows),
          'Table 8. Reference fragmentation: the same reads against a complete, a '
          'fragmented, and a scaffolded reference',
          'The 100-contig reference holds 43,707 of the complete genome\'s 43,735 '
          'anchors, so the genomic coordinate is the only variable. qc_pass is the '
          'share of estimates the fusion QC accepts -- note that it is HIGHEST '
          'where the estimates are worst. stationary_control is the estimate for '
          'the RUN_OUT sample, whose true log2(PTR) is ~0. Blank r or slope: the '
          'method returned nothing, or a constant.')


# --- Table 9: landmark-source comparison, one row per property ---------------
def table9():
    """Discussion landmark table from outline.md, with the F1/F2/F3/F4/F5/C5
    evidence filled in and unmeasured rows marked. Every row carries its
    evidence source."""
    rows = [
        ('landmark rule',
         'hash-prefix subsample of all 31-mers',
         'Type IIB recognition motifs, enumerated per reference',
         'Table 1; manifest.json'),
        ('coordinate enters estimator?',
         'no — sorted ranks',
         'yes — windowed V-fit',
         'Fig 6 (2x2 attribution)'),
        ('reference determinism',
         'deterministic under the same hash-prefix rule',
         'deterministic by motif enumeration',
         'both modes; the determinism claim was dropped'),
        ('internal-consistency strata',
         'hash-prefix split possible, but shares every systematic bias of the sketch',
         '~15 motif strata with heterogeneous biases (motif / GC / methylatable bases)',
         "Cochran's Q; Bsp24I-CjePI containment"),
        ('accuracy at matched density, >=1x',
         'identical',
         'identical',
         'F1 (density-matched sweep, this study)'),
        ('robustness at 0.5x',
         'loses 6/17 cells at the no-gradient gate (single stratum, unrescuable)',
         '15–16/17 cells via multi-strata fusion redundancy',
         'F1 (mechanism cells: per-window fill equal; fusion redundancy)'),
        ('usable GC range',
         'tied with panel; boundary GC <=30% at 0.5x',
         'same boundary (instrument limit, both modes)',
         'F2 (18 genomes, GC 25.4-72.0%)'),
        ('divergence geometry',
         'smooth: one SNP destroys ~26% of a read\'s 31-mers',
         'blocky: one SNP kills one site, neighbours intact',
         'retention at 0.1% substitutions: 94.6% vs 89.5% (4 enzymes); F4 rel-arm'),
        ('mismatch tolerance',
         'must lock mm=0: mm>=1 opens a rescue channel crediting ~1% (mm1) / ~2% (mm2) '
         'of observations to the wrong coordinate, independent of sequencing error',
         'mm<=2 safe: exact-first suppression + motif gate keep misassignment <=1.4e-4',
         'F3 (near-neighbor census + read-level test)'),
        ('wrong / absent reference',
         'silent — the gate refuses output',
         'loud — the estimate is emitted and wrong',
         'C5: reported_fraction 1.00 (denominator artifact) vs silent gate; '
         'common-denominator recall 3.8-13.8% vs 5.0-11.1%'),
        ('fragmented reference',
         'rank regression barely affected (r 0.889->0.827)',
         'gradient destroyed, QC-blind, scaffold-repairable',
         'Fig 9; F4: collapse identical across landmark types '
         '(slope 0.26-0.38, RMSE x40, winner\'s-curse residue)'),
        ('scaffolding',
         'works — placement is landmark-source-agnostic',
         'works — self-arm ties the sketch arm exactly',
         'F4: 98-100/100 contigs placed, order Spearman 1.0, both modes; '
         'sketch placement is a 1:1 port of scaffold.rs, not shipped'),
        ('density control',
         'scale parameter (content-agnostic)',
         'panel composition (biology-aware, per-clade tunable)',
         'Table 6 ranking; per-clade panel survey: not measured'),
        ('index cost at matched density',
         '~39 B/landmark if engineered into the same index format; merged FMH index '
         '(--mode fracminhash) not implemented',
         '~38.7-38.9 B/anchor; merged 300-genome DB builds 10M anchors in 86 s',
         'F5'),
        ('half-density option',
         'scale 200: ~124 GB GTDB projection, <=0.007 r loss at >=1x',
         'no equivalent knob (panel size is the knob; Table 7)',
         'F5'),
        ('wet-lab realization',
         'computational only',
         'physically produced by the 2bRAD protocol',
         'route B'),
        ('methylation-bias effect on yield',
         'n/a',
         'not measured',
         '—'),
    ]
    t = pd.DataFrame(rows, columns=['property', 'fracminhash_pilea',
                                    'enzyme_anchors_ours', 'evidence'])
    write('table9_landmark_comparison', t,
          'Table 9. Landmark-source comparison: property by mode, every row with '
          'its evidence',
          'Rows marked "not measured" have no data in this study and must not be '
          'read as null results. The mismatch-tolerance row is asymmetric for a '
          'measured reason (F3): the enzyme mode\'s exact-first counting rule '
          'suppresses the near-neighbor search that becomes the sketch mode\'s '
          'rescue channel at scale ~100. Sketch-side scaffolding and sketch '
          'mismatch tolerance rest on a rule-based port/simulation respectively, '
          'not on shipped CLI modes (F3/F4 REVIEW caveats).')


# --- Table 10: metagenome results, three datasets ----------------------------
def table10():
    blocks = []
    sun = DATA / 'sun_three_arm'
    c5 = DATA / 'c5_review'

    # (a) Sun cohort, cross-method agreement on WGS arms.
    f = sun / 'p2_review_agreement.tsv'
    if f.exists():
        d = pd.read_csv(f, sep='\t')
        d = d.rename(columns={'r': 'pearson_r', 'range_x': 'sk2bgrow_wgs_range',
                              'range_y': 'comparator_range'})
        order = {'A_default_vs_B': 0, 'A_gatesoff_vs_B': 1, 'B_vs_C': 2}
        d['ord'] = d['comparison'].map(order)
        a = (d[d['comparison'] != 'B_vs_C'].sort_values(['ord', 'sample'])
             .drop(columns='ord').reset_index(drop=True))
        blocks.append(('(a) Sun fecal cohort (PRJNA689204): cross-method agreement '
                       'on the WGS arms', a))

    # (b) Real 2bRAD arm: raw vs deduped against the WGS arm. Deduped rows are
    # the B_vs_C rows of the agreement tsv; raw rows exist only in REVIEW §6
    # (no tsv), and pass fewer species through the same gates.
    if f.exists():
        ded = (d[d['comparison'] == 'B_vs_C'].drop(columns='ord')
               .reset_index(drop=True))
        ded.insert(0, 'arm', 'real 2bRAD, deduped')
        raw = pd.DataFrame([
            ('real 2bRAD, raw', 'S01', 8, 0.449, 0.440, -0.23),
            ('real 2bRAD, raw', 'S06', 10, 0.797, 0.728, 0.60),
            ('real 2bRAD, raw', 'S07', 16, 0.734, 0.566, 0.62),
        ], columns=['arm', 'sample', 'n', 'pearson_r', 'ccc', 'bias'])
        b = pd.concat([ded, raw], ignore_index=True)[
            ['arm', 'sample', 'n', 'pearson_r', 'ccc', 'bias',
             'loa_lo', 'loa_hi', 'sk2bgrow_wgs_range', 'comparator_range']]
        blocks.append(('(b) Same cohort, in-silico vs real 2bRAD library (same '
                       'method, two preparations)', b))

    # (b, cont.) per-anchor capture efficiency, from sigma_icc.tsv.
    sf = sun / 'sigma_icc.tsv'
    if sf.exists():
        s = pd.read_csv(sf, sep='\t')
        eff = pd.DataFrame([
            ('sigma_eff, raw per-anchor dispersion', 1.53,
             'first per-anchor capture-efficiency measurement for 2bRAD'),
            ('ICC-correctable fraction, median (IQR)',
             f"{s['correctable_frac'].median():.3f} "
             f"({s['correctable_frac'].quantile(0.25):.3f}-"
             f"{s['correctable_frac'].quantile(0.75):.3f})",
             'between/(between+within), Poisson-corrected; 21 species with '
             '>=30 shared anchors across the 3 samples'),
            ('residual sigma_eff after e_prior, median',
             f"{s['resid_sigma'].median():.3f}",
             'within-batch; window averaging (/sqrt(99)) dilutes it to '
             '~0.07-0.10 on the gradient scale'),
        ], columns=['metric', 'value', 'note'])
        blocks.append(('(b, cont.) Per-anchor capture efficiency (route-B '
                       'characterization)', eff))

    # (c) C5 RBC application: recall three ways, QC validity, cost.
    rf = c5 / 'c5_recall_three_ways.tsv'
    cf = c5 / 'c5_common_denominator.tsv'
    qf = c5 / 'c5_coverage_control.tsv'
    if rf.exists():
        r = pd.read_csv(rf, sep='\t')
        sk = r[r['arm'] == 'sk2bgrow']
        pil = r[r['arm'] == 'pilea_default']
        n_common_qc = ''
        if cf.exists():
            c = pd.read_csv(cf, sep='\t')
            n_common_qc = (f"; sk2bgrow QC passes {int(c['n_common_qc_sk2b'].min())}-"
                           f"{int(c['n_common_qc_sk2b'].max())} of "
                           f"{int(c['n_common'].min())}-{int(c['n_common'].max())} on it")
        rho_frag = rho_cov = None
        if qf.exists():
            q = pd.read_csv(qf, sep='\t').set_index('feature')
            rho_frag = q.loc['n_contigs', 'rho_raw']
            rho_cov = q.loc['n_contigs', 'rho_partial_cov']
        c5t = pd.DataFrame([
            ('sk2bgrow reported_fraction',
             f"{sk['reported_fraction'].min():.2f} ({int(sk['n_rows'].iloc[0])}/"
             f"{int(sk['n_rows'].iloc[0])})",
             'no output gate; a denominator artifact, never a performance claim'),
            ('sk2bgrow estimate_fraction',
             f"{sk['estimate_fraction'].min():.3f}-{sk['estimate_fraction'].max():.3f}",
             '19 rows carry no PTR estimate (coverage NaN, 0 QC pass)'),
            ('sk2bgrow qc_recall',
             f"{100*sk['qc_recall'].min():.1f}-{100*sk['qc_recall'].max():.1f}% "
             f"({int(sk['n_qc_pass'].min())}-{int(sk['n_qc_pass'].max())} of 522)",
             'the same-denominator comparator to Pilea default'),
            ('Pilea default reported',
             f"{100*pil['reported_fraction'].min():.1f}-"
             f"{100*pil['reported_fraction'].max():.1f}% "
             f"({int(pil['n_rows'].min())}-{int(pil['n_rows'].max())} of 522)",
             'overlaps sk2bgrow qc_recall almost exactly'),
            ('common denominator',
             'Pilea-default MAGs are a strict subset of sk2bgrow outputs'
             + n_common_qc,
             'c5_common_denominator.tsv'),
            ('QC pass vs MAG fragmentation',
             (f"Spearman rho {rho_frag:.2f} vs n_contigs; {rho_cov:.2f} controlling "
              f"mean coverage" if rho_frag is not None else
              'see c5_coverage_control.tsv'),
             'QC works on real data; coverage is the strongest single predictor '
             'but the fragmentation effect survives the control'),
            ('cost vs Pilea, per sample',
             '89.5-240.8x (measured, c5_cost_per_sample.tsv) -> ~10-25x after --max-mismatch 1',
             'mm=1: lookup 63x faster (collapsed middle seed), 5.3% anchors lost, '
             '0 genomes lost; C5 REVIEW §3/§4f'),
        ], columns=['metric', 'value', 'note'])
        blocks.append(('(c) C5 rotating biological contactor (PRJNA974210, 9 '
                       'samples, 522 MAGs)', c5t))

    note = ('Agreement coefficients are over species x sample units. In (a) the '
            'comparator is Pilea at its shipped gates; the observed ranges span '
            'only ~1.4 log2, so the low CCC (0.25-0.37) is dynamic-range '
            'deflation, not broken concordance; the gates-off rows (r ~0) are '
            'the control showing Pilea\'s gate does real work. In (b) exact '
            'read-level deduplication before counting flattens the PTR dynamic '
            'range (deduped r 0.30-0.60 vs raw 0.45-0.80); the cap-dedup '
            'sensitivity sweep found no robust intermediate (cap=2 is best in '
            'S07, worst in S06), so the process recommendation is to not '
            'exactly dedup 2bRAD libraries before counting; a GLM on deduped '
            'counts does not recover the signal (r -0.03-0.41). sigma_eff is '
            'within-batch (3 libraries, one study/one centre; cross-lab '
            'transfer awaits the Hou cohort) and measured on the single enzyme '
            '(BcgI) present in the Sun libraries. In (c) the 1.00 recall '
            'headline is reported_fraction (no gate), not accuracy.')
    write_blocks('table10_metagenome', blocks,
                 'Table 10. Metagenome results: Sun fecal-cohort concordance, the '
                 'real-2bRAD dedup finding, and the C5 RBC application', note)


if __name__ == '__main__':
    print('generating tables ->', OUT)
    table1(); table2(); table3(); table5(); table6(); table7(); table8()
    table9(); table10()
    print('done')
