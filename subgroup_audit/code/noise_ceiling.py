"""Classical measurement-noise ceiling, with cluster-bootstrap uncertainty.

This is an approximation requiring unbiased additive errors and independent
replicate components. It is NOT a measured maximum for sign accuracy.
"""
import argparse,hashlib
from pathlib import Path
import pandas as pd,numpy as np
from common import ROOT,plan,write_json
from audit_stats import weights,noise,pearson
from audit_engine import group_sum

def audit(d,draws=10000):
 rows=[]
 for dataset,g in d.loc[d.primary].groupby('dataset'):
  q=g.loc[np.isfinite(g.noise_variance)&g.noise_variance.ge(0)&np.isfinite(g.effect)].copy().reset_index(drop=True)
  if len(q)<20:
   rows.append({'dataset':dataset,'eligible_variants':len(g),'replicate_noise_variants':len(q),'mean_effect_reliability':np.nan,'Pearson_noise_ceiling':np.nan,'ceiling_ci95_low':np.nan,'ceiling_ci95_high':np.nan,'status':'NOT IDENTIFIABLE: full benchmark lacks per-variant replicate effect estimates' if dataset=='Kircher' else 'NOT IDENTIFIABLE: too few replicate estimates','approximation':'No SE manufactured from publication p-values'});continue
  q['cluster']=q.genomic_block;clusters=sorted(q.cluster.unique());y=q.effect.to_numpy(float);v=q.noise_variance.to_numpy(float)
  s=group_sum(q,{'n':np.ones(len(q)),'y':y,'yy':y*y,'noise':v},clusters);r,c=noise(s.sum(axis=0));seed=int(hashlib.sha256(('913071noise'+dataset).encode()).hexdigest()[:8],16);w=weights(len(clusters),seed,draws);rb,cb=noise(w@s);valid=np.isfinite(cb);ci=np.quantile(cb[valid],[.025,.975]) if valid.any() else [np.nan,np.nan]
  rec={'dataset':dataset,'eligible_variants':len(g),'replicate_noise_variants':len(q),'clusters':len(clusters),'mean_effect_reliability':float(r),'Pearson_noise_ceiling':float(c),'ceiling_ci95_low':ci[0],'ceiling_ci95_high':ci[1],'bootstrap_identifiable_fraction':float(valid.mean()),'reliability_ci95_low':np.nanquantile(rb,.025),'reliability_ci95_high':np.nanquantile(rb,.975),'status':'IDENTIFIABLE APPROXIMATION' if np.isfinite(c) and valid.mean()>=.95 else 'UNIDENTIFIED OR UNSTABLE','approximation':'Abell: 3 paired transfections; Ulirsch: RNA/DNA crossed components, unmodeled library covariance remains'}
  for model in ['enformer','borzoi']:
   col=model+'_DNASE_mean'
   if col not in q:continue
   x=q[col].to_numpy(float);ps=group_sum(q,{'n':np.ones(len(q)),'x':x,'y':y,'xx':x*x,'yy':y*y,'xy':x*y},clusters)
   raw=float(pearson(ps.sum(axis=0)));corr=raw/c if np.isfinite(c) else np.nan;cr=pearson(w@ps)/cb;ok=np.isfinite(cr);ci2=np.quantile(cr[ok],[.025,.975]) if ok.any() else [np.nan,np.nan]
   rec[model+'_observed_Pearson']=raw;rec[model+'_noise_corrected_Pearson']=corr;rec[model+'_corrected_ci95_low']=ci2[0];rec[model+'_corrected_ci95_high']=ci2[1];rec[model+'_correction_status']='UNSTABLE |r|>1, not clipped' if abs(corr)>1 else 'Conditional additive-error approximation' if np.isfinite(corr) else 'UNIDENTIFIED'
  rows.append(rec)
 return pd.DataFrame(rows)

def library_agreement(d):
 from scipy.stats import pearsonr,spearmanr
 rows=[]
 for stem in ['LDLR','SORT1']:
  a=d.loc[d.dataset.eq('Kircher')&d.condition.eq(stem)&d.assay_eligible,['variant_id','effect']]
  b=d.loc[d.dataset.eq('Kircher')&d.condition.eq(stem+'.2')&d.assay_eligible,['variant_id','effect']]
  q=a.merge(b,on='variant_id',suffixes=('_1','_2'),validate='one_to_one')
  rows.append({'element':stem,'paired_library_variants':len(q),'Pearson_library_agreement':pearsonr(q.effect_1,q.effect_2).statistic if len(q)>3 else np.nan,'Spearman_library_agreement':spearmanr(q.effect_1,q.effect_2).statistic if len(q)>3 else np.nan,'note':'Two independent library measurements for this element only; not a full Kircher noise ceiling and no independent-variant p-value claim'})
 return pd.DataFrame(rows)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--variants',required=True);p.add_argument('--out',default=str(ROOT/'results/public_data_noise_ceiling.csv'));a=p.parse_args();plan();d=pd.read_csv(a.variants);z=audit(d);z.to_csv(a.out,index=False);print(z.to_string(index=False))
