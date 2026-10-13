"""Subgroup metadata from frozen coordinates and actual selected track IDs."""
import json,re
from pathlib import Path
import pandas as pd,numpy as np
from common import ROOT,MODELS,sha

def expected_cell(r):
 if bool(r.cross_species_proxy):return 'mouse assay / human model'
 if r.dataset=='Ulirsch':return 'K562'
 if r.dataset=='Abell':return 'GM12878'
 if r.condition in ['TERT-GBM','TERT-GAa','TERT-GSc']:return 'SF7996'
 s=str(r.cell)
 if s.startswith('LNCaP'):return 'LNCaP'
 if s.startswith('HEK293T'):return 'HEK293T'
 return s

def annotate(d,base):
 aliases=json.loads((ROOT/'config/cell_aliases.json').read_text());manifest=json.loads((ROOT/'annotations/annotation_manifest.json').read_text());anno=ROOT/'annotations/variant_ccre_annotations.csv.gz'
 if sha(anno)!=manifest['subset_sha256']:raise ValueError('cCRE subset checksum changed')
 a=pd.read_csv(anno)
 if not a.variant_id.is_unique:raise ValueError('Duplicate cCRE variant annotations')
 if not set(d.variant_id)<=set(a.variant_id):raise ValueError('Input universe differs; do not silently label new variants unannotated')
 d=d.merge(a[['variant_id','regulatory_type','ccre_signatures','ccre_ids']],on='variant_id',validate='many_to_one');rows=[]
 for model in MODELS:
  meta=pd.read_csv(Path(base)/'data/model_metadata'/(model+'_targets.tsv'),sep='\t')
  if 'Unnamed: 0' in meta:meta=meta.rename(columns={'Unnamed: 0':'index'})
  sel=[]
  for p in sorted((Path(base)/'scores'/model).glob('tracks*.csv')):
   z=pd.read_csv(p);sel.extend(z.loc[z.modality.eq('DNASE')].to_dict('records'))
  if not sel:raise ValueError('Actual DNASE track selections missing for '+model)
  mapping={}
  for r in sel:
   key=(r['element'],r['condition']);ids=tuple(int(float(x)) for x in str(r['track_indices']).split(';') if x and x!='nan')
   if key in mapping and mapping[key]!=ids:raise ValueError('Track selections differ across shards')
   mapping[key]=ids
  cellstatus={}
  for _,r in d[['dataset','element','condition','cell','cross_species_proxy']].drop_duplicates().iterrows():
   key=(r.element,r.condition);ids=mapping.get(key)
   if ids is None:raise ValueError('No actual selected track IDs: '+str(key))
   m=meta.set_index('index').loc[list(ids)];text=m.description.fillna('').astype(str);cell=expected_cell(r);pattern=aliases.get(cell)
   if cell=='mouse assay / human model':status='unmatched'
   elif not len(text) or text.eq('').any():status='unknown'
   elif pattern is None:status='unknown'
   elif text.str.contains(pattern,case=False,regex=True).all():status='matched'
   elif cell=='HEL 92.1.7' and text.str.contains(r'(?<![A-Za-z0-9])HEL(?![A-Za-z0-9])',case=False,regex=True).all():status='unknown'
   else:status='unmatched'
   cellstatus[(r.dataset,r.element,r.condition)]=status
   rows.append({'model':model,'dataset':r.dataset,'element':r.element,'condition':r.condition,'expected_cell':cell,'cell_match':status,'selected_DNASE_track_indices':';'.join(map(str,ids)),'selected_descriptions':' || '.join(text)})
  d[model+'_cell_match']=[cellstatus[(r.dataset,r.element,r.condition)] for r in d.itertuples()]
 return d,pd.DataFrame(rows)
