"""Reproducible helpers. Original Wahab artifacts are read-only inputs."""
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parents[1]
MODELS=['enformer','borzoi'];MODS=['DNASE','ATAC','CAGE','H3K27ac']
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 return h.hexdigest()
def clean(x):
 import numpy as np
 if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [clean(v) for v in x]
 if isinstance(x,np.generic):x=x.item()
 if isinstance(x,float) and not np.isfinite(x):return None
 return x
def write_json(x,p):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix('.partial');tmp.write_text(json.dumps(clean(x),indent=2,allow_nan=False)+'\n');tmp.replace(p)
def plan():
 p=ROOT/'config/preregistered_plan.json';lock=json.loads((ROOT/'config/plan_lock.json').read_text())
 if sha(p)!=lock['sha256']:raise ValueError('Locked subgroup plan changed. Use a separate new version, never edit this run.')
 return json.loads(p.read_text())
def input_reader(base):
 """Validate complete native-model journals without rewriting N-padding hashes.

 Code identities are retained as provenance; identical data, protocol and
 checkpoint model specs are required. No predictor is executed or imported.
 """
 import pandas as pd,numpy as np
 base=Path(base).resolve();d=pd.read_csv(base/'data/variants.csv');p=base/'config/protocol.json'
 if not d.row_id.is_unique:raise ValueError('Duplicate source row IDs')
 source=sha(base/'data/variants.csv');protocol=sha(p);provenance=[]
 for model in MODELS:
  rows=[];specs=set()
  for s in range(8):
   folder=base/'scores'/model;j=folder/f'shard{s:02d}.jsonl';state=folder/f'shard{s:02d}.json'
   z=json.loads(state.read_text());i=z['identity']
   if not z.get('complete') or not z.get('real_model_scores'):raise ValueError('Incomplete native model shard '+model+str(s))
   if i['source_rows_sha256']!=source or i['protocol_sha256']!=protocol:raise ValueError('Source/protocol identity mismatch')
   if i['model']!=model or i['shards']!=8 or i['shard']!=s:raise ValueError('Shard identity mismatch')
   specs.add(json.dumps(i['model_spec'],sort_keys=True));provenance.append({'model':model,'shard':s,'state_sha256':sha(state),'journal_sha256':sha(j),'identity':i,'recorded_device':z.get('device')})
   for line in j.read_text().splitlines():
    r=json.loads(line)
    if r.get('model')!=model or not r.get('real_model_scores'):raise ValueError('Non-native or wrong model record')
    rows.append(r)
  if len(specs)!=1:raise ValueError('Mixed weight/context model specs')
  f=pd.DataFrame(rows)
  if not f.row_id.is_unique or set(f.row_id)!=set(d.row_id):raise ValueError('Missing/duplicate model rows')
  for c in ['DNASE_mean','DNASE_std']:
   if not np.isfinite(f[c]).all():raise ValueError('Primary DNASE scores are missing')
  columns=['row_id']+[m+'_'+s for m in MODS for s in ['mean','std','n']]
  f=f[columns].rename(columns={c:model+'_'+c for c in columns if c!='row_id'})
  d=d.merge(f,on='row_id',validate='one_to_one')
 # The original primary audit, not an optional model-dependent analysis.
 summary=json.loads((base/'results/audit_summary.json').read_text())
 if summary.get('AlphaGenome_in_primary_fit') or not summary.get('real_model_scores'):raise ValueError('Original primary fit is not the required open-model analysis')
 e=pd.read_csv(base/'results/external_predictions.csv')
 hashes=base/'results/result_CSV_sha256.json'
 if hashes.exists():
  hh=json.loads(hashes.read_text())
  if hh.get('external_predictions.csv')!=sha(base/'results/external_predictions.csv'):raise ValueError('Original external prediction export changed')
 pool=d.loc[d.primary & d.dataset.isin(['Abell','Ulirsch'])]
 if not e.row_id.is_unique or set(e.row_id)!=set(pool.row_id):raise ValueError('Original policy export differs from validation universe')
 a=pool[['row_id','enformer_DNASE_mean']].merge(e[['row_id','forecast']],on='row_id',validate='one_to_one')
 np.testing.assert_allclose(a.enformer_DNASE_mean,a.forecast,rtol=0,atol=1e-12)
 if not np.isfinite(e.risk_raw).all():raise ValueError('Missing original policy risks')
 risk=e.set_index('row_id').risk_raw
 d['original_enformer_risk']=d.row_id.map(risk)
 return d,{'base':str(base),'source_sha256':source,'protocol_sha256':protocol,'original_primary_summary':summary,'shards':provenance,'known_implementation_fix':'User reports N-padding correction at chromosome boundaries; recorded code hashes preserved, no scoring changes made by this package.'}
