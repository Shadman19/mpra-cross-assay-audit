"""Reconstruct Kircher mapped contrasts directly from the checked public table.

Raw GP1BA is displayed as GP1BB per its promoter coordinates. HBG2 was not a
separate element. Multiple conditions remain acquired; source training uses one
fixed canonical condition per element. HBB/HBG1 and the MYC regions share groups.
"""
import numpy as np
import pandas as pd
from common import ROOT
CANONICAL=['BCL11A','F9','FOXE1','GP1BA','HBB','HBG1','HNF4A','IRF4','IRF6','LDLR','MSMB','MYCrs11986220','MYCrs6983267','PKLR-48h','RET','SORT1','TCF7L2','TERT-HEK','UC88','ZFAND3','ZRSh-13']
CELLS={'BCL11A':'HEL 92.1.7','F9':'HepG2','FOXE1':'HeLa','GP1BB':'HEL 92.1.7','HBB':'HEL 92.1.7','HBG1':'HEL 92.1.7','HNF4A':'HEK293T','IRF4':'SK-MEL-28','IRF6':'HaCaT','LDLR':'HepG2','MSMB':'HEK293T','MYCrs11986220':'LNCaP + DHT','MYCrs6983267':'HEK293T + LiCl','PKLR':'K562','RET':'Neuro-2a (mouse cells, human insert)','SORT1':'HepG2','TCF7L2':'MIN6 (mouse cells, human insert)','TERT':'HEK293T','UC88':'Neuro-2a (mouse cells, human insert)','ZFAND3':'MIN6 (mouse cells, human insert)','ZRS':'NIH/3T3 + HOXD13 (mouse cells, human insert)'}
def family(s):
    for stem in ['TERT','SORT1','LDLR','PKLR']:
        if s.startswith(stem):return stem
    if s.startswith('ZRSh'):return 'ZRS'
    return 'GP1BB' if s=='GP1BA' else s

def reconstruct(root=ROOT):
    d=pd.read_csv(root/'data/raw/kircher_elements.tsv.gz',sep='\t',dtype={'Chrom':str},low_memory=False)
    d=d.loc[d.Release.eq('GRCh38')].copy()
    mapped=d.Ref.isin(list('ACGT')) & d.Alt.isin(list('ACGT')+['-']) & d.Ref.ne(d.Alt) & d.Chrom.notna() & d.Pos.notna()
    un=d.loc[~mapped].copy();un['reason']='No source coordinates/canonical alleles'
    un.to_csv(root/'data/kircher_unmapped_rows.csv',index=False)
    d=d.loc[mapped].copy();d['element']=d.Element.map(family)
    d['holdout_group']=d.element.replace({'HBB':'GLOBIN','HBG1':'GLOBIN','MYCrs11986220':'MYC','MYCrs6983267':'MYC'})
    d['assay_eligible']=d.Barcodes.ge(10)&np.isfinite(d.Coefficient)
    d['primary']=d.Element.isin(CANONICAL)&d.Alt.ne('-')&d.assay_eligible
    d['chrom']='chr'+d.Chrom;d['pos']=d.Pos.astype(int);d['ref']=d.Ref;d['alt']=d.Alt;d['effect']=d.Coefficient.astype(float)
    d['variant_id']=d.chrom+':'+d.pos.astype(str)+':'+d.ref+'>'+d.alt;d['row_id']=d.Element+'|'+d.variant_id
    d['assay_cell']=d.element.map(CELLS);d['cross_species_proxy']=d.element.isin(['RET','UC88','TCF7L2','ZFAND3','ZRS'])
    d.loc[d.Element.isin(['TERT-GBM','TERT-GAa','TERT-GSc']),'assay_cell']='Glioblastoma SF7996, with stated siRNA condition'
    cols=['row_id','variant_id','Element','element','holdout_group','chrom','pos','ref','alt','effect','Barcodes','DNA','RNA','pValue','assay_eligible','primary','assay_cell','cross_species_proxy']
    assert d.row_id.is_unique
    d[cols].to_csv(root/'data/kircher_prepared.csv',index=False)
    return d[cols]
