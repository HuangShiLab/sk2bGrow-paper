#!/usr/bin/env python3
from __future__ import annotations
import argparse, math, re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score, average_precision_score, balanced_accuracy_score, f1_score
from sklearn.model_selection import LeaveOneOut, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

def numeric(df, cols):
    for c in cols:
        if c in df: df[c] = pd.to_numeric(df[c], errors='coerce')
    return df

def label_for(sample: str, dataset: str):
    s = sample.upper()
    if dataset == 'PRJNA946904':
        if s.startswith('MIBC'): return 'MIBC'
        if s.startswith('NMIBC'): return 'NMIBC'
    if dataset == 'PRJNA1005427':
        if sample.startswith('C-'): return 'cancer'
        if sample.startswith('B-'): return 'benign'
    if dataset == 'CRA019236':
        if sample.startswith('MST.'): return 'maternal_feces'
        if sample.startswith('MEC.'): return 'infant_meconium'
        if sample.startswith('MK.T2.'): return 'breast_milk_T2'
        if sample.startswith('MK.T4.'): return 'breast_milk_T4'
    return None

def read_outputs(root: Path):
    parts=[]
    for acc in ['PRJNA946904','PRJNA1005427','CRA019236']:
        for p in sorted((root/'per_sample_sp'/acc).glob('*/noGC_glm/output.tsv')):
            sample=p.parent.parent.name
            d=pd.read_csv(p,sep='\t')
            d['sample_id']=sample; d['dataset']=acc
            parts.append(d)
    if not parts:
        raise SystemExit('no output.tsv files found')
    d=pd.concat(parts,ignore_index=True)
    cols=['PTR','log2(PTR)','se','coverage','dispersion','fraction','containment','n_anchors','n_windows','ori_confidence','enzyme_fit_rate']
    d=numeric(d,cols)
    d['log2_ptr']=d['log2(PTR)']
    d['label']=d.apply(lambda r: label_for(r.sample_id,r.dataset),axis=1)
    return d

def robust_mask(d):
    return d['log2_ptr'].notna() & (d['n_windows']>=5) & (d['n_anchors']>=50) & (d['ori_confidence']>=0.5) & (d['se']<=1) & d['log2_ptr'].between(0,3)

def wide(d, value_cols, suffix):
    z=d.set_index(['dataset','sample_id','taxonomy'])[value_cols].unstack('taxonomy')
    z.columns=[f'{a}__{b}{suffix}' for a,b in z.columns]
    z.index=z.index.get_level_values('sample_id')
    return z

def select_train_columns(X,y,min_nonmissing=5):
    keep=[]
    for c in X.columns:
        if X.loc[y.index.difference([y.name] if False else [])].empty: continue
        if X[c].notna().sum() >= min_nonmissing and X[c].nunique(dropna=True)>1:
            keep.append(c)
    fallback=list(X.columns[X.columns.str.endswith('__missing')])
    return keep or fallback or [X.columns[0]]

def cv_arm(X,y,seed=20260926,n_trees=300):
    # Imputation uses median and adds explicit missing indicators. No zero imputation for PTR.
    pipe=Pipeline([
      ('impute',SimpleImputer(strategy='median',add_indicator=True)),
      ('scale',StandardScaler()),
      ('rf',RandomForestClassifier(n_estimators=n_trees,n_jobs=-1,class_weight='balanced',random_state=seed)),
    ])
    splitter=LeaveOneOut() if len(y)<=50 else StratifiedKFold(5,shuffle=True,random_state=seed)
    prob=[]; pred=[]; test_samples=[]
    classes=np.array(sorted(y.unique()))
    for train,test in splitter.split(np.zeros(len(y)),y):
        ytr,yte=y.iloc[train],y.iloc[test]
        Xtr=X.iloc[train]
        keep=select_train_columns(Xtr,ytr,min_nonmissing=max(3,min(8,len(ytr)//3)))
        pipe.fit(Xtr.loc[:,keep],ytr)
        p=pipe.predict_proba(X.iloc[test].loc[:,keep])
        prob.append(p); pred.extend(pipe.predict(X.iloc[test].loc[:,keep]))
        test_samples.extend(y.index[test])
    prob=np.vstack(prob); pred=np.array(pred); yy=(y==classes[-1]).astype(int) if len(classes)==2 else pd.get_dummies(y).reindex(columns=classes,fill_value=0).to_numpy()
    m={'n_samples':len(y),'n_features_before':X.shape[1]}
    if len(classes)==2:
        out=pd.DataFrame({'sample':test_samples,'y_true':y.to_numpy(),'prob_positive':prob[:,1],'pred':pred})
        m.update(auroc=roc_auc_score(yy,prob[:,1]),auprc=average_precision_score(yy,prob[:,1]))
        return m,out
    else:
        m.update(auroc=roc_auc_score(pd.get_dummies(y).reindex(columns=classes,fill_value=0),prob,multi_class='ovr',average='macro'))
        m.update(auprc=np.nan)
    m.update(balanced_accuracy=balanced_accuracy_score(y,pred),macro_f1=f1_score(y,pred,average='macro'))
    out=pd.DataFrame({'sample':test_samples,'y_true':y.to_numpy(),'pred':pred})
    out=pd.concat([out,pd.DataFrame(prob,columns=[f'prob_{c}' for c in classes])],axis=1)
    return m,out

def classify(d,outdir):
    metrics=[]; predictions=[]
    for ds,g in d.groupby('dataset'):
        samples=g[['sample_id','label']].drop_duplicates().dropna()
        samples=samples[samples.label.notna()]
        if len(samples)<6: continue
        y=samples.set_index('sample_id')['label']
        gidx=g.set_index('sample_id')
        growth=gidx[gidx['log2_ptr'].notna()][['log2_ptr','se']]
        support=gidx[gidx['coverage'].notna()][['coverage','n_windows','n_anchors','containment','ori_confidence']]
        A=wide(gidx.reset_index(),['coverage','n_windows','n_anchors','containment','ori_confidence'],'__support')
        # Only robust PTR estimates enter the growth feature matrix. Remaining taxa are missing,
        # and are represented by explicit missingness indicators rather than zero imputation.
        B=wide(gidx.loc[robust_mask(gidx)].reset_index(),['log2_ptr','se'],'__ptr')
        A=A.reindex(y.index); B=B.reindex(y.index)
        # Explicit missingness indicators; median imputation is used later, never zero imputation.
        def add_missing_indicators(X):
            miss=X.isna()
            miss=miss.loc[:,miss.any()]
            miss.columns=[c+'__missing' for c in miss.columns]
            return pd.concat([X,miss.astype(float)],axis=1) if not miss.empty else X
        A=add_missing_indicators(A); B=add_missing_indicators(B)
        arms={'A_taxonomic_support':A,'B_growth_PTR':B,'C_combined':pd.concat([A,B],axis=1)}
        for name,X in arms.items():
            m,out=cv_arm(X,y); m.update(dataset=ds,arm=name,n_labels=y.value_counts().to_dict())
            metrics.append(m)
            out['dataset']=ds; out['arm']=name; predictions.append(out)
        # Save matrices for reproducibility.
        A.to_csv(outdir/f'{ds}_support_matrix.tsv.gz',sep='\t',compression='gzip')
        B.to_csv(outdir/f'{ds}_ptr_matrix.tsv.gz',sep='\t',compression='gzip')
    return pd.DataFrame(metrics),pd.concat(predictions,ignore_index=True) if predictions else pd.DataFrame()

def bootstrap_auc_ci(y,s,n_boot=5000,seed=20260926):
    y=np.asarray(y); s=np.asarray(s); rng=np.random.default_rng(seed); vals=[]
    for _ in range(n_boot):
        idx=rng.integers(0,len(y),len(y)); yy=y[idx]
        if len(np.unique(yy))<2: continue
        vals.append(roc_auc_score(yy,s[idx]))
    return (float(np.quantile(vals,0.025)),float(np.quantile(vals,0.975))) if vals else (np.nan,np.nan)

def paired_auc_tests(pred_df,y,n_perm=5000,seed=20260926):
    rows=[]
    datasets=pred_df.dataset.unique()
    for ds in datasets:
        z=pred_df[pred_df.dataset==ds]
        classes=sorted(z.y_true.unique())
        if len(classes)!=2: continue
        yy=(z.groupby('sample').y_true.first()==classes[-1]).astype(int)
        probs={arm:z[z.arm==arm].set_index('sample').prob_positive for arm in z.arm.unique()}
        rng=np.random.default_rng(seed)
        for contrast,(a,b) in {'B-A':('A_taxonomic_support','B_growth_PTR'),'C-A':('A_taxonomic_support','C_combined')}.items():
            if a not in probs or b not in probs: continue
            pa,pb=probs[a].reindex(yy.index),probs[b].reindex(yy.index)
            delta=roc_auc_score(yy,pb)-roc_auc_score(yy,pa); null=[]
            for _ in range(n_perm):
                perm=pd.Series(rng.permutation(yy.to_numpy()),index=yy.index)
                null.append(roc_auc_score(perm,pb)-roc_auc_score(perm,pa))
            p=(np.sum(np.abs(np.asarray(null))>=abs(delta))+1)/(n_perm+1)
            rows.append({'dataset':ds,'contrast':contrast,'auc_delta':delta,'p_permutation':p})
    return pd.DataFrame(rows)

def bootstrap_metrics(pred_df,n_boot=5000,seed=20260926):
    rows=[]
    for (ds,arm),z in pred_df.groupby(['dataset','arm']):
        yy=z.set_index('sample').y_true; classes=sorted(yy.unique())
        if len(classes)!=2: continue
        s=z.set_index('sample').prob_positive
        lo,hi=bootstrap_auc_ci((yy==classes[-1]).astype(int),s,n_boot,seed)
        rows.append({'dataset':ds,'arm':arm,'auc_ci_lower':lo,'auc_ci_upper':hi,'positive_class':classes[-1]})
    return pd.DataFrame(rows)

def bh(p):
    p=np.asarray(p); n=len(p); o=np.argsort(p); ranked=p[o]*(n/np.arange(1,n+1)); q=np.minimum.accumulate(ranked[::-1])[::-1]; z=np.empty(n); z[o]=np.minimum(q,1); return z

def contrasts(d,outdir):
    rows=[]
    def fit(ds,name,ga,gb,gcol='label'):
        g=d[(d.dataset==ds)&d[gcol].notna()&d.log2_ptr.notna()&robust_mask(d)]
        for tax,z in g.groupby('taxonomy'):
            a=z.loc[z[gcol]==ga,'log2_ptr'].dropna(); b=z.loc[z[gcol]==gb,'log2_ptr'].dropna()
            if len(a)>=3 and len(b)>=3 and (a.nunique()>1 or b.nunique()>1):
                try: stat,p=mannwhitneyu(b,a,alternative='two-sided')
                except ValueError: continue
                rows.append({'dataset':ds,'contrast':f'{gb}_vs_{ga}','taxonomy':tax,'n_'+ga:len(a),'n_'+gb:len(b),'median_'+ga:a.median(),'median_'+gb:b.median(),'median_log2ptr_effect':b.median()-a.median(),'p':p})
    fit('PRJNA946904','bladder','MIBC','NMIBC')
    fit('PRJNA1005427','ovarian','benign','cancer')
    fit('CRA019236','maternal_vs_meconium','maternal_feces','infant_meconium')
    fit('CRA019236','milk_T4_vs_T2','breast_milk_T2','breast_milk_T4')
    r=pd.DataFrame(rows)
    if not r.empty:
        for _,z in r.groupby('contrast'): r.loc[z.index,'q']=bh(z['p'].to_numpy())
    return r.sort_values(['contrast','q','p'])

def paired_milk(d,outdir):
    z=d[(d.dataset=='CRA019236')&d.label.isin(['breast_milk_T2','breast_milk_T4'])&robust_mask(d)].copy()
    z['subject']=z.sample_id.str.replace(r'^MK\.T[24]\.','',regex=True)
    rows=[]
    for tax,g in z.groupby('taxonomy'):
        w=g.pivot(index='subject',columns='label',values='log2_ptr').dropna()
        if len(w)>=5 and (w['breast_milk_T4']-w['breast_milk_T2']).abs().sum()>0:
            stat,p=wilcoxon(w['breast_milk_T4'],w['breast_milk_T2'])
            rows.append({'taxonomy':tax,'n_pairs':len(w),'median_T2':w['breast_milk_T2'].median(),'median_T4':w['breast_milk_T4'].median(),'median_effect_T4_T2':(w['breast_milk_T4']-w['breast_milk_T2']).median(),'p':p})
    r=pd.DataFrame(rows)
    if r.empty:
        return pd.DataFrame(columns=['taxonomy','n_pairs','median_T2','median_T4','median_effect_T4_T2','p','q'])
    if not r.empty: r['q']=bh(r.p.to_numpy())
    return r.sort_values('q')

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',required=True,type=Path)
    ap.add_argument('--outdir',required=True,type=Path)
    a=ap.parse_args(); a.outdir.mkdir(parents=True,exist_ok=True)
    d=read_outputs(a.root)
    d['robust']=robust_mask(d)
    d.to_csv(a.outdir/'all_output_long.tsv.gz',sep='\t',index=False,compression='gzip')
    # Yield summary
    y=d.groupby(['dataset','sample_id']).agg(n_species_tested=('taxonomy','size'),n_finite_ptr=('log2_ptr',lambda x:x.notna().sum()),n_robust_ptr=('robust','sum'),median_n_windows=('n_windows','median')).reset_index()
    y.to_csv(a.outdir/'per_sample_feature_yield.tsv',sep='\t',index=False)
    ysum=y.groupby('dataset').agg(n_samples=('sample_id','size'),mean_finite_ptr=('n_finite_ptr','mean'),median_finite_ptr=('n_finite_ptr','median'),mean_robust_ptr=('n_robust_ptr','mean'),median_robust_ptr=('n_robust_ptr','median')).reset_index()
    ysum.to_csv(a.outdir/'dataset_feature_yield_summary.tsv',sep='\t',index=False)
    m,pred_df=classify(d,a.outdir); m.to_csv(a.outdir/'classification_loocv_metrics.tsv',sep='\t',index=False)
    if not pred_df.empty:
        pred_df.to_csv(a.outdir/'classification_loocv_predictions.tsv',sep='\t',index=False)
        paired_auc_tests(pred_df,pd.Series(dtype=str)).to_csv(a.outdir/'paired_auc_tests.tsv',sep='\t',index=False)
        bootstrap_metrics(pred_df).to_csv(a.outdir/'auc_bootstrap_ci.tsv',sep='\t',index=False)
    c=contrasts(d,a.outdir); c.to_csv(a.outdir/'group_contrasts_log2ptr.tsv',sep='\t',index=False)
    p=paired_milk(d,a.outdir); p.to_csv(a.outdir/'breast_milk_T4_vs_T2_paired.tsv',sep='\t',index=False)
    print(ysum.to_string(index=False)); print('\nCLASSIFICATION\n',m.to_string(index=False)); print('\nTOP CONTRASTS\n',c.head(30).to_string(index=False))
if __name__=='__main__': main()
