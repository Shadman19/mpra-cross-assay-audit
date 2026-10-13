"""Descriptive assay-view coverage, without extra tests or selector fitting."""
from pathlib import Path
import pandas as pd,numpy as np
from common import ROOT

def coverage(d):
 rows=[]
 for dataset,g in d.groupby('dataset'):
  screens=g.effect.to_numpy()[:,None] if dataset=='Kircher' else g[[f'audit_screen_{v}' for v in 'ABC']].to_numpy(float)
  valid=np.isfinite(screens);strong=valid&(abs(screens)>=1);weak=valid&(abs(screens)>=.5)&(abs(screens)<1);negligible=valid&(abs(screens)<.5)
  rows.append({'dataset':dataset,'variants':len(g),'strong_views':int(strong.sum()),'weak_active_views':int(weak.sum()),'negligible_noisy_views':int(negligible.sum()),'variants_with_strong_view':int(strong.any(axis=1).sum()),'variants_with_weak_active_view':int(weak.any(axis=1).sum()),'variants_with_negligible_view':int(negligible.any(axis=1).sum()),'variants_in_multiple_bins':int(((strong.any(axis=1).astype(int)+weak.any(axis=1).astype(int)+negligible.any(axis=1).astype(int))>1).sum()),'note':'Bins can overlap across replicate views; negligible views not added to sign-target tests. No new p-values or label-based selection.'})
 return pd.DataFrame(rows)
if __name__=='__main__':
 d=pd.read_csv(ROOT/'results/annotated_audit_variants.csv.gz');coverage(d).to_csv(ROOT/'results/assay_effect_view_coverage.csv',index=False)
