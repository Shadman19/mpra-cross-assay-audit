"""Prepare real allele contrasts and replicate-aware audit targets.

Ulirsch: centered inserts, pooled barcodes, 6 RNA libraries and 2 DNA libraries.
Abell: Base single-nucleotide contrasts in 1KG pools, three paired RNA/DNA reps.
No model output is available or consulted at this stage.
"""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
import py2bit
from pyliftover import LiftOver
from common import ROOT,verify_sources,load_protocol,atomic_json
COMP=str.maketrans('ACGT','TGCA')
def paired_effect(rref,ralt,dref,dalt):
    return np.log2((ralt+.5)/(rref+.5))-np.log2((dalt+.5)/(dref+.5))
def pad_ids(chr_,pos,ref,alt):return str(chr_)+':'+str(int(pos))+':'+str(ref)+'>'+str(alt)
def prepare(root=ROOT,reference=None):
    protocol=load_protocol(root);verify_sources(root)
    refpath=Path(reference) if reference else root/'runtime/hg38.2bit'
    if not refpath.exists():raise FileNotFoundError('Run code/download_reference.py first, or supply --reference /path/hg38.2bit')
    genome=py2bit.open(str(refpath));lo=LiftOver(str(root/'data/raw/hg19ToHg38.over.chain.gz'))
    records=[];reps=[];exclusions=[]
    from prepare_kircher import reconstruct
    k=reconstruct(root)
    for r in k.itertuples():
        rid='Kircher|'+r.row_id
        records.append({'row_id':rid,'legacy_row_id':r.row_id,'dataset':'Kircher','condition':r.Element,'element':r.element,'fit_group':r.holdout_group,'chrom':r.chrom,'pos':int(r.pos),'ref':r.ref,'alt':r.alt,'variant_id':r.variant_id,'cell':r.assay_cell,'primary':bool(r.primary),'effect':r.effect,'noise_variance':np.nan,'audit_screen_A':np.nan,'audit_eval_A':np.nan,'audit_screen_B':np.nan,'audit_eval_B':np.nan,'audit_screen_C':np.nan,'audit_eval_C':np.nan,'assay_eligible':bool(r.assay_eligible),'source_replication':'No variant-level replicate estimates for most elements; never manufacture SEs','cross_species_proxy':bool(r.cross_species_proxy)})
    # Ulirsch barcode rows are repeated observations of the same allele/insert.
    # Centered window is selected in advance, not whichever window looks best.
    u=pd.read_csv(root/'data/raw/ulirsch_counts.tsv.gz',sep='\t',low_memory=False)
    u=u.loc[u.clean.eq('var') & u.bot.eq('1/2') & u.top.eq('1/2')].copy()
    keys=['chr','pos','ref','alt','type'];countcols=[c for c in u if c.startswith('K562_')]
    u=u.groupby(keys,dropna=False)[countcols].sum().reset_index()
    for key,g in u.groupby(['chr','pos','ref','alt'],sort=True):
        chrom0,pos0,ref0,alt0=key;chrom0='chr'+str(chrom0).removesuffix('.0');pos0=int(pos0)
        if ref0 not in 'ACGT' or alt0 not in 'ACGT' or ref0==alt0:continue
        if set(g.type)!=set(['Ref','Mut']):continue
        a=g.set_index('type');rr=a.loc['Ref',[f'K562_CTRL_minP_RNA{i}' for i in range(1,7)]].to_numpy(float);ra=a.loc['Mut',[f'K562_CTRL_minP_RNA{i}' for i in range(1,7)]].to_numpy(float)
        dr=a.loc['Ref',['K562_minP_DNA1','K562_minP_DNA2']].to_numpy(float);da=a.loc['Mut',['K562_minP_DNA1','K562_minP_DNA2']].to_numpy(float)
        origin=pad_ids(chrom0,pos0,ref0,alt0)
        mapped=lo.convert_coordinate(chrom0,pos0-1)
        if len(mapped)!=1:exclusions.append({'dataset':'Ulirsch','source':origin,'reason':'Ambiguous/unmapped hg19 to hg38 coordinate'});continue
        chrom,p0,strand,_=mapped[0];ref,alt=ref0,alt0
        if strand=='-':ref=ref.translate(COMP);alt=alt.translate(COMP)
        vid=pad_ids(chrom,p0+1,ref,alt);rid='Ulirsch|'+vid
        # Both RNA and DNA uncertainty are included. A shared DNA denominator
        # otherwise falsely makes RNA replicate agreement look like full truth.
        rna=np.log2((ra+.5)/(rr+.5));dna=np.log2((da+.5)/(dr+.5))
        effect=float(rna.mean()-dna.mean());noise=float(rna.var(ddof=1)/6+dna.var(ddof=1)/2)
        quality=bool(np.minimum(dr,da).min()>=10)
        records.append({'row_id':rid,'legacy_row_id':origin,'dataset':'Ulirsch','condition':'K562_CTRL','element':'Ulirsch_variant','fit_group':chrom+':'+str(int((p0+1)//1000000)),'chrom':chrom,'pos':int(p0+1),'ref':ref,'alt':alt,'variant_id':vid,'cell':'K562','primary':quality,'effect':effect,'noise_variance':noise,'audit_screen_A':float(rna[:3].mean()-dna[0]),'audit_eval_A':float(rna[3:].mean()-dna[1]),'audit_screen_B':float(rna[3:].mean()-dna[1]),'audit_eval_B':float(rna[:3].mean()-dna[0]),'audit_screen_C':np.nan,'audit_eval_C':np.nan,'assay_eligible':quality,'source_replication':'6 RNA / 2 DNA libraries; covariance-independence approximation','cross_species_proxy':False})
        for i,x in enumerate(rna):reps.append({'row_id':rid,'component':'RNA_logratio','replicate':i+1,'effect':float(x)})
        for i,x in enumerate(dna):reps.append({'row_id':rid,'component':'DNA_logratio','replicate':i+1,'effect':float(x)})
    # Abell contains haplotype perturbations and reused Tewhey controls. Only
    # Base SNVs in independent 1KG pools are primary, with every exclusion logged.
    a=pd.read_csv(root/'data/raw/abell_counts.tsv.gz',sep='\t',low_memory=False)
    conv=pd.read_csv(root/'data/raw/abell_hg_conversion.tsv',sep='\t').drop_duplicates()
    if conv.duplicated(['hg19_chrom','hg19_pos']).any():raise ValueError('Ambiguous Abell conversion table')
    a=a.merge(conv,left_on=['chrom','pos'],right_on=['hg19_chrom','hg19_pos'],how='left',validate='many_to_one')
    for r in a.itertuples():
        source=str(r.VarID)
        if r.Pool not in ['1KGTest','1KGControl'] or r.haploFlag!='Base' or len(str(r.ref))!=1 or len(str(r.alt))!=1 or r.ref not in 'ACGT' or r.alt not in 'ACGT':
            exclusions.append({'dataset':'Abell','source':source,'reason':'Non-Base haplotype, reused Sabeti pool, or non-SNV'});continue
        if pd.isna(r.hg38_pos):exclusions.append({'dataset':'Abell','source':source,'reason':'No author hg38 conversion'});continue
        chrom=str(r.hg38_chrom);pos=int(r.hg38_pos);vid=pad_ids(chrom,pos,r.ref,r.alt);rid='Abell|'+source+'|'+vid
        rr=np.array([getattr(r,f'cDNA{i}_ref') for i in range(1,4)],float);ra=np.array([getattr(r,f'cDNA{i}_alt') for i in range(1,4)],float)
        dr=np.array([getattr(r,f'plasmid{i}_ref') for i in range(1,4)],float);da=np.array([getattr(r,f'plasmid{i}_alt') for i in range(1,4)],float)
        e=paired_effect(rr,ra,dr,da);quality=bool(np.minimum(dr,da).min()>=10)
        rec={'row_id':rid,'legacy_row_id':source,'dataset':'Abell','condition':'GM12878','element':'Abell_variant','fit_group':chrom+':'+str(pos//1000000),'chrom':chrom,'pos':pos,'ref':r.ref,'alt':r.alt,'variant_id':vid,'cell':'GM12878','primary':quality,'effect':float(e.mean()),'noise_variance':float(e.var(ddof=1)/3),'assay_eligible':quality,'source_replication':'3 paired RNA/DNA transfections','cross_species_proxy':False}
        for j,label in enumerate('ABC'):
            rec['audit_screen_'+label]=float(e[j]);rec['audit_eval_'+label]=float(np.delete(e,j).mean())
        records.append(rec)
        for i,x in enumerate(e):reps.append({'row_id':rid,'component':'paired_logratio_effect','replicate':i+1,'effect':float(x)})
    d=pd.DataFrame(records)
    if not d.row_id.is_unique:raise ValueError('Duplicate source contrasts; do not silently average them')
    # Reference validation is done for every contrast. Unmatched library alleles
    # are acquired as explicit contrasts but excluded from external primary tests.
    matches=[]
    for r in d.itertuples():
        try:matches.append(genome.sequence(r.chrom,int(r.pos)-1,int(r.pos)-1+len(r.ref)).upper()==r.ref)
        except Exception:matches.append(False)
    d['reference_match']=matches
    d.loc[d.dataset.ne('Kircher') & ~d.reference_match,'primary']=False
    # Independent experiments are insufficient if the same locus reappears.
    # Lock a conservative 1 Mb block exclusion between Kircher and validation,
    # then Ulirsch and Abell; no outcome or model score chooses this exclusion.
    d['genomic_block']=d.chrom+':'+(d.pos//1000000).astype(str)
    kircher_blocks=set(d.loc[d.dataset.eq('Kircher'),'genomic_block'])
    ul_blocks=set(d.loc[d.dataset.eq('Ulirsch'),'genomic_block'])
    d['independent_locus']=~d.genomic_block.isin(kircher_blocks)
    d.loc[d.dataset.eq('Abell'),'independent_locus'] &= ~d.loc[d.dataset.eq('Abell'),'genomic_block'].isin(ul_blocks)
    d.loc[d.dataset.ne('Kircher') & ~d.independent_locus,'primary']=False
    # Several Abell Base constructs encode the identical genomic contrast.
    # Keep a source-ID-defined canonical QC-eligible contrast, never count
    # repeated constructs as additional independent variants or pick by effect.
    canonical=d.loc[d.dataset.eq('Abell') & d.primary].sort_values('row_id').drop_duplicates('variant_id').row_id
    d.loc[d.dataset.eq('Abell') & ~d.row_id.isin(canonical),'primary']=False
    d['tile_start']=(d.pos.sub(1)//2048*2048).astype(int)
    d['context_id']=d.chrom+':'+d.tile_start.astype(str)
    d.to_csv(root/'data/variants.csv',index=False)
    pd.DataFrame(reps).to_csv(root/'data/replicates.csv',index=False)
    pd.DataFrame(exclusions).to_csv(root/'data/preparation_exclusions.csv',index=False)
    qc=[]
    for dataset,g in d.groupby('dataset'):
        p=g.loc[g.primary]
        qc.append({'dataset':dataset,'acquisition_rows':len(g),'primary_variants':len(p),'unique_primary_variants':p.variant_id.nunique(),'primary_blocks':p.genomic_block.nunique(),'observed_abs_effect_ge_0_5':int(p.effect.abs().ge(.5).sum()),'reference_mismatches':int((~g.reference_match).sum()),'replicate_information':g.source_replication.iloc[0]})
    pd.DataFrame(qc).to_csv(root/'results/data_QC.csv',index=False)
    atomic_json({'actual_public_data_prepared':True,'model_predictions_generated':False,'model_results':'Pending','dataset_QC':qc},root/'results/preparation_status.json')
    genome.close();print(pd.DataFrame(qc).to_string(index=False));return d
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--reference');a=p.parse_args();prepare(reference=a.reference)
