"""Rebuild the bundled, cell-agnostic cCRE annotations on public coordinates.

One-off annotation preparation only. Statistical runs need neither pyBigWig
nor the50MB BigBed, because the coordinate-matched subset is bundled.
"""
import argparse,json
from pathlib import Path
import pandas as pd
from common import ROOT,sha,write_json,plan
URL='https://hgdownload.soe.ucsc.edu/gbdb/hg38/encode3/ccre/encodeCcreCombined.bb'
EXPECTED='42a9cf761ed87e2ce7b26d60a84b9b623848a2a8c2de9c84c246eb4c6d38417c'
def build(variants,bigbed,out=None):
 import pyBigWig,numpy as np
 plan()
 if sha(bigbed)!=EXPECTED:raise ValueError('ENCODE3 BigBed hash differs from frozen annotation release')
 d=pd.read_csv(variants,usecols=['variant_id','chrom','pos']).drop_duplicates('variant_id');bb=pyBigWig.open(str(bigbed));records=[]
 for chrom,g in d.groupby('chrom',sort=True):
  positions=np.sort(g.pos.unique().astype(int)-1);types={int(p):set() for p in positions};ids={int(p):set() for p in positions}
  if chrom in bb.chroms():
   start=max(0,int(positions.min()));end=min(bb.chroms(chrom),int(positions.max())+1)
   entries=bb.entries(chrom,start,end) if end>start else None
   for a,b,rest in entries or []:
    parts=rest.split('\t');labels=set(parts[6].split(',')) # BED9+6 ccre field.
    i=np.searchsorted(positions,a);j=np.searchsorted(positions,b)
    for p in positions[i:j]:types[int(p)].update(labels);ids[int(p)].add(parts[0])
  for r in g.itertuples():
   t=types[int(r.pos)-1];prom='PLS' in t;enh=bool(t & {'pELS','dELS'})
   kind='ambiguous' if prom and enh else 'promoter' if prom else 'enhancer' if enh else 'other' if t else 'unannotated'
   records.append({'variant_id':r.variant_id,'chrom':chrom,'pos':int(r.pos),'regulatory_type':kind,'ccre_signatures':';'.join(sorted(t)),'ccre_ids':';'.join(sorted(ids[int(r.pos)-1]))})
 bb.close();target=Path(out) if out else ROOT/'annotations/variant_ccre_annotations.csv.gz'
 pd.DataFrame(records).sort_values('variant_id').to_csv(target,index=False,compression={'method':'gzip','mtime':0})
 write_json({'source_url':URL,'source_sha256':EXPECTED,'source_bytes':49809776,'anonymous_download_verified':True,'annotation':'ENCODE3 hg38/UCSC encodeCcreCombined, cell-agnostic BED9+6; NOT an experimental MPRA construct-type label','coordinate_source_sha256':sha(variants),'subset_sha256':sha(target),'subset_variants':len(records),'coordinate_convention':'1-based variant POS converted to0-based half-open overlap','unknown':'Other/unannotated/ambiguous never relabeled as enhancer'},ROOT/'annotations/annotation_manifest.json')
 print('Bundled public-coordinate cCRE annotation:',len(records),'variants')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--variants',required=True);p.add_argument('--bigbed',required=True);a=p.parse_args();build(a.variants,a.bigbed)
