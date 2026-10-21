#!/usr/bin/env python3
from __future__ import annotations
import argparse, re
from pathlib import Path
import pandas as pd
from scipy.stats import pearsonr, spearmanr

NUM=['PTR','log2_ptr','se','coverage','dispersion','fraction','containment','n_anchors','n_windows','ori_confidence','enzyme_fit_rate']
TIME={'9':'9AM','11':'11AM','13':'1PM','17':'5PM'}
def num(df,cols):
    for c in cols:
        if c in df: df[c]=pd.to_numeric(df[c],errors='coerce')
    return df

def robust_mask(d):
    return (d.log2_ptr.notna() & (d.n_windows>=5) & (d.n_anchors>=50)
            & (d.ori_confidence>=0.5) & (d.se<=1) & d.log2_ptr.between(0,3))

def strict(d):
    return d.pass_qc.astype(str).str.lower().isin(['true','1'])

def key_from_wms(x): return re.sub(r'^saliva_WMS__','',x)
def key_from_2brad(library): return re.sub(r'^2bRAD_','',library)

def add_demo(d,keycol='sample_key'):
    d['host']=d[keycol].str.extract(r'^(S\d+)_',expand=False)
    d['timepoint_code']=d[keycol].str.extract(r'^S\d+_(.+)$',expand=False)
    d['timepoint']=d['timepoint_code'].map(TIME).fillna(d['timepoint_code'])
    return d

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',required=True,type=Path)
    ap.add_argument('--routeB-long',required=True,type=Path)
    ap.add_argument('--routeB-metadata',required=True,type=Path)
    ap.add_argument('--outdir',required=True,type=Path)
    a=ap.parse_args(); a.outdir.mkdir(parents=True,exist_ok=True)
    # WMS route-A rows
    wparts=[]
    for method in ('auto','glm'):
        for p in sorted((a.root/'per_sample').glob(f'*/noGC_{method}/output.tsv')):
            d=pd.read_csv(p,sep='\t').rename(columns={'log2(PTR)':'log2_ptr'}); d=num(d,NUM)
            sid=p.parent.parent.name
            d=d.assign(sample_id=sid,sample_key=key_from_wms(sid),assay='WMS',
                       wms_method=method,run_method=d.get('method'))
            wparts.append(d)
    if not wparts: raise SystemExit('no WMS outputs')
    w=pd.concat(wparts,ignore_index=True)
    w=add_demo(w)
    w['robust']=robust_mask(w); w['pass_qc_bool']=strict(w)
    # routeB 2bRAD rows, saliva subset only
    rb=pd.read_csv(a.routeB_long,sep='\t'); rb=num(rb,NUM)
    md=pd.read_csv(a.routeB_metadata,sep='\t')
    md=md[md.library_name.str.match(r'^2bRAD_S(?:9|11|13|17)_',na=False)].copy()
    md['sample_key']=md.library_name.map(key_from_2brad)
    keys=set(md.sample_key)
    rb=rb[(rb.project=='PRJNA1131785') & rb.sample_id.isin(md.sample_id)].copy()
    rb=rb.merge(md[['sample_id','sample_key','library_name']],on='sample_id',how='left',validate='many_to_one')
    rb=rb[rb.sample_key.isin(keys)].copy()
    rb=add_demo(rb)
    rb=rb.assign(assay='2bRAD',wms_method='glm',run_method='glm')
    rb['robust']=rb.robust
    rb['pass_qc_bool']=strict(rb)
    cols=['sample_id','sample_key','host','timepoint_code','timepoint','assay','wms_method','genome','taxonomy',
          'PTR','log2_ptr','se','coverage','dispersion','fraction','containment','n_anchors','n_windows','ori_confidence',
          'enzyme_fit_rate','pass_qc_bool','robust','run_method']
    long=pd.concat([w[cols],rb[cols]],ignore_index=True)
    long.to_csv(a.outdir/'combined_2brad_wms_long.tsv.gz',sep='\t',index=False,compression='gzip')
    # WMS yield
    wy=w.assign(robust=robust_mask(w),strict=strict(w)).groupby(['sample_key','host','timepoint','wms_method'],as_index=False).agg(
        n_species=('genome','size'),n_finite=('log2_ptr','size'),n_robust=('robust','size'),n_strict=('strict','size'),
        median_coverage=('coverage','median'),median_n_anchors=('n_anchors','median'))
    wy.to_csv(a.outdir/'wms_per_sample_yield.tsv',sep='\t',index=False)
    # routeB yield for matched 32
    ry=rb.groupby(['sample_key','host','timepoint'],as_index=False).agg(
        n_species=('genome','size'),n_finite=('log2_ptr','size'),n_robust=('robust','size'),n_strict=('pass_qc','size'))
    ry.to_csv(a.outdir/'matched_2brad_per_sample_yield.tsv',sep='\t',index=False)
    # Wide species-sample pairing
    r=rb[['sample_key','genome','taxonomy','host','timepoint_code','timepoint','log2_ptr','PTR','se','coverage','n_windows','n_anchors','ori_confidence','robust','pass_qc']].copy()
    r=r.rename(columns={'log2_ptr':'log2ptr_2bRAD','PTR':'PTR_2bRAD','se':'se_2bRAD','coverage':'coverage_2bRAD','n_windows':'nwin_2bRAD','n_anchors':'nanchors_2bRAD','ori_confidence':'ori_2bRAD','robust':'robust_2bRAD','pass_qc':'strict_2bRAD'})
    for method in ('auto','glm'):
        z=w[w.wms_method==method][['sample_key','genome','taxonomy','host','timepoint_code','timepoint','log2_ptr','PTR','se','coverage','n_windows',
        'n_anchors','ori_confidence','pass_qc']].copy()
        z=z.assign(robust=robust_mask(z),strict=strict(z))
        z=z.rename(columns={'log2_ptr':f'log2ptr_WMS_{method}','PTR':f'PTR_WMS_{method}','se':f'se_WMS_{method}','coverage':f'coverage_WMS_{method}','n_windows':f'nwin_WMS_{method}','n_anchors':f'nanchors_WMS_{method}','ori_confidence':f'ori_WMS_{method}','robust':f'robust_WMS_{method}','strict':f'strict_WMS_{method}'})
        r=r.merge(z.drop(columns=['host','timepoint_code','timepoint']),on=['sample_key','genome','taxonomy'],how='outer',suffixes=('','_tmp'))
    for method in ('auto','glm'):
        lp2='log2ptr_2bRAD'; lpw=f'log2ptr_WMS_{method}'
        r[f'delta_log2ptr_WMS_minus_2bRAD_{method}']=r[lpw]-r[lp2]
    r.to_csv(a.outdir/'paired_species_sample_wide.tsv.gz',sep='\t',index=False,compression='gzip')
    # correlations
    rows=[]
    for method in ('auto','glm'):
        lpw=f'log2ptr_WMS_{method}'; rw=f'robust_WMS_{method}'; sw=f'strict_WMS_{method}'
        sub=r[[lpw,'log2ptr_2bRAD','robust_2bRAD','strict_2bRAD',rw,sw]].copy()
        specs=[('all_finite',sub.log2ptr_2bRAD.notna()&sub[lpw].notna()),
               ('both_exploratory_robust',sub.robust_2bRAD.fillna(False)&sub[rw].fillna(False)),
               ('2bRAD_robust_only',sub.robust_2bRAD.fillna(False)),
               ('WMS_robust_only',sub[rw].fillna(False)),
               ('both_strict_qc',sub.strict_2bRAD.fillna(False)&sub[sw].fillna(False))]
        for name,mask in specs:
            z=sub.loc[mask.fillna(False),[lpw,'log2ptr_2bRAD']].dropna()
            rec={'wms_method':method,'filter':name,'n_pairs':len(z)}
            if len(z)>=3 and z.iloc[:,0].nunique()>1 and z.iloc[:,1].nunique()>1:
                rec['pearson_r']=pearsonr(z.iloc[:,1],z.iloc[:,0]).statistic
                rec['spearman_rho']=spearmanr(z.iloc[:,1],z.iloc[:,0]).statistic
                rec['median_abs_delta']=float((z.iloc[:,0]-z.iloc[:,1]).abs().median())
            rows.append(rec)
    pd.DataFrame(rows).to_csv(a.outdir/'correlation_by_filter.tsv',sep='\t',index=False)
    print('WMS rows',len(w),'2bRAD rows',len(rb),'paired wide',len(r))
    print(pd.DataFrame(rows).to_string(index=False))
if __name__=='__main__': main()
