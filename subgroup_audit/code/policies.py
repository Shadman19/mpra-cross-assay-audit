"""Open-model-only CPU risk fitting. No foundation model or API is invoked.

Feature ordering reproduces the original v2 Enformer recipe. Borzoi uses the
same source-only recipe with model identities swapped. Hyperparameters frozen.
"""
import hashlib
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.isotonic import IsotonicRegression
SEED=20261012

def features(d,model,scale=None):
 other='borzoi' if model=='enformer' else 'enformer';f=d[model+'_DNASE_mean'].to_numpy(float);o=d[other+'_DNASE_mean'].to_numpy(float)
 if not np.isfinite(f).all() or not np.isfinite(o).all():raise ValueError('Missing primary forecasts')
 if scale is None:scale={'self':max(1e-6,float(np.std(f))),'other':max(1e-6,float(np.std(o)))}
 keys=[d.dataset,d.genomic_block]
 fs=d[model+'_DNASE_mean'].groupby(keys).transform(lambda v:float(np.std(v,ddof=0))).to_numpy()/scale['self']
 os=d[other+'_DNASE_mean'].groupby(keys).transform(lambda v:float(np.std(v,ddof=0))).to_numpy()/scale['other']
 cols=[np.log1p(abs(f)/scale['self']),d[model+'_DNASE_std'].to_numpy(float)/scale['self'],abs(f/scale['self']-o/scale['other']),fs,os]
 for name in [model,other]:
  for mod in ['ATAC','CAGE','H3K27ac']:
   v=d[name+'_'+mod+'_mean'].to_numpy(float);cols.extend([np.nan_to_num(v,nan=0),np.isfinite(v).astype(float)])
 return np.column_stack(cols),scale

def balanced(g):
 counts=pd.Series(g).value_counts();w=np.array([1/counts[v] for v in g]);return w/w.mean()

def fit(source,model):
 groups=sorted(source.fit_group.unique(),key=lambda v:hashlib.sha256((str(SEED)+v).encode()).hexdigest());n=int(np.ceil(.25*len(groups)))
 calgroups=groups[:n];fitgroups=groups[n:];f=source.loc[source.fit_group.isin(fitgroups)];c=source.loc[source.fit_group.isin(calgroups)]
 x,scale=features(f,model);mask=f.effect.abs().ge(.5).to_numpy();forecast=f[model+'_DNASE_mean'].to_numpy();y=(np.sign(forecast[mask])!=np.sign(f.effect.to_numpy()[mask])).astype(int)
 if mask.sum()<20:raise ValueError('Source fit has fewer than20 effectful observations')
 est=HistGradientBoostingClassifier(max_iter=120,max_leaf_nodes=7,min_samples_leaf=30,l2_regularization=5,learning_rate=.07,random_state=SEED)
 constant=float(y.mean()) if np.unique(y).size<2 else None
 if constant is None:est.fit(x[mask],y,sample_weight=balanced(f.fit_group.to_numpy()[mask]))
 iso=None;cx,_=features(c,model,scale);raw=np.full(len(c),constant) if constant is not None else est.predict_proba(cx)[:,1];cm=c.effect.abs().ge(.5).to_numpy()
 cy=(np.sign(c[model+'_DNASE_mean'].to_numpy()[cm])!=np.sign(c.effect.to_numpy()[cm])).astype(int)
 if cm.sum() and np.unique(raw[cm]).size>1:iso=IsotonicRegression(out_of_bounds='clip').fit(raw[cm],cy,sample_weight=balanced(c.fit_group.to_numpy()[cm]))
 return {'estimator':est,'constant':constant,'calibrator':iso,'scale':scale,'model':model,'fit_groups':fitgroups,'calibration_groups':calgroups,'features':'17 fixed output-only features; no measured validation values','source_fit_effectful':int(mask.sum()),'source_calibration_effectful':int(cm.sum())}

def predict(d,p):
 x,_=features(d,p['model'],p['scale']);raw=np.full(len(d),p['constant']) if p['constant'] is not None else p['estimator'].predict_proba(x)[:,1]
 return raw

def source_lofo(source,model):
 result=pd.Series(np.nan,index=source.index);manifest=[]
 for group in sorted(source.fit_group.unique()):
  train=source.loc[source.fit_group.ne(group)];held=source.loc[source.fit_group.eq(group)];p=fit(train,model)
  if group in p['fit_groups'] or group in p['calibration_groups']:raise AssertionError('Held source family leaked')
  result.loc[held.index]=predict(held,p);manifest.append({'heldout':group,'fit_groups':p['fit_groups'],'calibration_groups':p['calibration_groups']})
 return result,manifest
