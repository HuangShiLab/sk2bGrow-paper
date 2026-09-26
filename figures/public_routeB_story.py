#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.metrics import roc_curve, auc
from audit_panel_alignment import require_matplotlib_panel_alignment

root=Path('analysis/public_routeB_external/results')
outdir=Path('figures/public_routeB_story')
outdir.mkdir(parents=True,exist_ok=True)
mpl.rcParams.update({
    'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans','sans-serif'],
    'pdf.fonttype':42,'svg.fonttype':'none','font.size':7,'axes.titlesize':7,
    'axes.labelsize':7,'xtick.labelsize':6,'ytick.labelsize':6,
    'axes.spines.right':False,'axes.spines.top':False,'axes.linewidth':0.7,
    'legend.frameon':False
})
colors={'A':'#4C78A8','B':'#E45756','C':'#72B7B2'}
arm_names={'A':'Support','B':'PTR','C':'Combined'}
ds_names={'PRJNA1005427':'Ovary','PRJNA946904':'Bladder','CRA019236':'Maternal'}
yield_df=pd.read_csv(root/'dataset_feature_yield_summary.tsv',sep='\t')
per_sample=pd.read_csv(root/'per_sample_feature_yield.tsv',sep='\t')
metrics=pd.read_csv(root/'classification_loocv_metrics.tsv',sep='\t')
boot=pd.read_csv(root/'auc_bootstrap_ci.tsv',sep='\t')
merged=metrics.merge(boot,on=['dataset','arm'])
pred=pd.read_csv(root/'classification_loocv_predictions.tsv',sep='\t')

fig,axs=plt.subplots(1,3,figsize=(7.087,2.35),constrained_layout=True)
# A robust yield
order=['PRJNA946904','PRJNA1005427','CRA019236']
labels=[ds_names[x] for x in order]
data=[per_sample.loc[per_sample.dataset==x,'n_robust_ptr'].to_numpy() for x in order]
bp=axs[0].boxplot(data,tick_labels=labels,patch_artist=True,widths=.55,medianprops={'color':'#333333','linewidth':1.0},
                  boxprops={'facecolor':'#E8EEF7','linewidth':.7},whiskerprops={'linewidth':.7},capprops={'linewidth':.7},flierprops={'markersize':2})
rng=np.random.default_rng(7)
for i,x in enumerate(data):
    axs[0].scatter(rng.uniform(i+1-.12,i+1+.12,len(x)),x,s=6,color='#2F5F8A',alpha=.65,linewidth=0)
axs[0].set_ylabel('Robust PTR estimates per sample')
axs[0].set_ylim(bottom=0)
axs[0].text(-.16,1.03,'a',transform=axs[0].transAxes,fontsize=9,fontweight='bold',va='bottom')

# B AUC + CI
pos=[]
for i,ds in enumerate(order):
    for j,arm in enumerate(['A_taxonomic_support','B_growth_PTR','C_combined']):
        r=merged[(merged.dataset==ds)&(merged.arm==arm)].iloc[0]
        x=i+(j-1)*.18
        color=list(colors.values())[j]
        yerr=np.array([[r.auroc-max(0,r.auc_ci_lower)],[min(1,r.auc_ci_upper)-r.auroc]])
        axs[1].errorbar(x,r.auroc,yerr=yerr,fmt='o',ms=3,capsize=1.6,elinewidth=.8,color=color)
        pos.append((x,r.auroc,color))
axs[1].set_xticks(range(3)); axs[1].set_xticklabels(labels)
axs[1].set_ylabel('LOOCV AUROC'); axs[1].set_ylim(.35,1.02)
handles=[Line2D([],[],marker='o',ls='',color=colors[k],label=arm_names[k],ms=3.2) for k in colors]
axs[1].legend(handles=handles,loc='lower left',ncol=1,handletextpad=.3,columnspacing=.7,borderaxespad=0)
axs[1].text(-.16,1.03,'b',transform=axs[1].transAxes,fontsize=9,fontweight='bold',va='bottom')

# C ROC curves, ovarian
roc_rows=[]
ax=axs[2]
for arm,label,color in [('A_taxonomic_support','Support',colors['A']),('B_growth_PTR','PTR',colors['B'])]:
    z=pred[(pred.dataset=='PRJNA1005427')&(pred.arm==arm)].groupby('sample').first()
    yy=(z.y_true=='cancer').astype(int); score=z.prob_positive
    fpr,tpr,_=roc_curve(yy,score); a=auc(fpr,tpr)
    ax.plot(fpr,tpr,color=color,lw=1.1,label=f'{label} (AUC {a:.2f})')
    roc_rows.append(pd.DataFrame({'fpr':fpr,'tpr':tpr,'arm':label}))
ax.plot([0,1],[0,1],ls='--',lw=.6,color='#777777')
ax.set_xlabel('False-positive rate'); ax.set_ylabel('True-positive rate')
ax.set_title('Ovarian: cancer vs benign',pad=2)
ax.legend(loc='lower right',borderaxespad=0,handlelength=1.2)
ax.text(-.16,1.03,'c',transform=ax.transAxes,fontsize=9,fontweight='bold',va='bottom')
# Source data
yield_df.to_csv(outdir/'panelA_dataset_summary.csv',index=False)
per_sample.to_csv(outdir/'panelA_per_sample.csv',index=False)
merged.to_csv(outdir/'panelB_auc.csv',index=False)
pd.concat(roc_rows,ignore_index=True).to_csv(outdir/'panelC_roc.csv',index=False)
require_matplotlib_panel_alignment(fig,json_out=str(outdir/'figure.alignment.json'),overlay_svg=str(outdir/'figure.alignment.svg'),tolerance_pt=1.5,gutter_tolerance_pt=1.5,strict=True)
fig.savefig(outdir/'public_routeB_story.svg',bbox_inches='tight')
fig.savefig(outdir/'public_routeB_story.pdf',bbox_inches='tight')
fig.savefig(outdir/'public_routeB_story.tiff',dpi=600,bbox_inches='tight')
