"""Publication-quality estimates, explicit uncertainty and unestimable strata."""
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
COLORS={'enformer':'#086E76','borzoi':'#BA563D'}
def build_figures(tests,noise,points,out):
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':320})
 def save(f,name):
  f.tight_layout()
  for suffix in ['png','pdf','svg']:f.savefig(out/(name+'.'+suffix),bbox_inches='tight')
  plt.close(f)
 q=tests.loc[tests.formal_external & tests.analysis_mode.eq('global') & tests.test_kind.eq('gain') & tests.dimension.ne('all')].copy()
 for ds,h in q.groupby('dataset'):
  f,a=plt.subplots(figsize=(8.2,6));h=h.sort_values(['dimension','level','model']).reset_index(drop=True)
  for i,r in h.iterrows():
   if not r.evaluable or not np.isfinite(r.estimate):
    a.text(.02,i,'Not estimable',transform=a.get_yaxis_transform(),color='gray',va='center',fontsize=8);continue
   x=r.estimate*100;lo=r.ci95_low*100;hi=r.ci95_high*100
   a.errorbar(x,i,xerr=[[max(0,x-lo)],[max(0,hi-x)]],fmt='o',mfc=COLORS[r.model] if r.p_Holm<.05 else 'white',mec=COLORS[r.model],ecolor=COLORS[r.model],capsize=3)
  a.set_yticks(np.arange(len(h)),h.dimension.str.replace('_',' ')+' / '+h.level+' / '+h.model);a.axvline(0,ls='--',color='gray');a.set_xlabel('Error reduction versus magnitude (percentage points)');a.set_title(ds+': original global selection, then stratification\n95% pointwise CIs; filled symbols pass70-test Holm correction',fontsize=10);save(f,'figure_subgroups_'+ds)
 f,a=plt.subplots(figsize=(7,4))
 for i,ds in enumerate(['Abell','Ulirsch']):
  for m,shift in [('enformer',-.12),('borzoi',.12)]:
   z=points.loc[points.dataset.eq(ds)&points.model.eq(m)].iloc[0];x=z.Pearson;a.errorbar(i+shift,x,yerr=[[max(0,x-z.ci95_low)],[max(0,z.ci95_high-x)]],fmt='o',color=COLORS[m],label=m if i==0 else None,capsize=4)
  n=noise.loc[noise.dataset.eq(ds)].iloc[0]
  if np.isfinite(n.Pearson_noise_ceiling):a.plot([i-.28,i+.28],[n.Pearson_noise_ceiling]*2,'--',color='gray',label='Approximate noise ceiling' if i==0 else None)
 a.set_xticks([0,1],['Abell','Ulirsch']);a.set(ylabel='Pearson correlation',title='Observed point fidelity and replicate-derived ceiling');a.legend(fontsize=8);save(f,'figure_noise_ceiling')
 f,a=plt.subplots(figsize=(6.8,4));q=tests.loc[tests.formal_external&tests.test_kind.eq('interaction')&tests.analysis_mode.eq('global')].copy().reset_index(drop=True)
 for i,r in q.iterrows():
  if not r.evaluable:continue
  x=r.estimate*100;a.errorbar(x,i,xerr=[[max(0,(r.estimate-r.ci95_low)*100)],[max(0,(r.ci95_high-r.estimate)*100)]],fmt='o',color=COLORS[r.model],capsize=3)
 a.set_yticks(np.arange(len(q)),q.dataset+' / '+q.model+' / '+q.dimension);a.axvline(0,ls='--',color='gray');a.set_xlabel('Difference in filtering gain (percentage points)');a.set_title('Prespecified subgroup interactions; missing contrasts unestimable');save(f,'figure_interactions')
