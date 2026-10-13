"""Locked subgroup contrasts. Measured-effect groups never drive selection."""
import hashlib
import numpy as np
import pandas as pd
from common import MODELS
from audit_stats import weights,gain,pearson,infer
LEVELS={'regulatory_type':['promoter','enhancer'],'effect_strength':['strong','weak'],'cell_match':['matched','unmatched']}

def errors(d,model,strength='all',source=False):
 forecast=d[model+'_DNASE_mean'].to_numpy(float)
 if source:screens=d.effect.to_numpy()[:,None];truth=screens.copy()
 else:screens=d[[f'audit_screen_{v}' for v in 'ABC']].to_numpy(float);truth=d[[f'audit_eval_{v}' for v in 'ABC']].to_numpy(float)
 mask=np.isfinite(screens)&np.isfinite(truth)&(truth!=0)&(abs(screens)>=.5)
 if strength=='strong':mask &= abs(screens)>=1
 if strength=='weak':mask &= abs(screens)<1
 bad=(np.sign(forecast[:,None])!=np.sign(truth)).astype(float);n=mask.sum(axis=1)
 return np.divide((bad*mask).sum(axis=1),n,out=np.full(len(d),np.nan),where=n>0)

def select(d,score,mask=None):
 """Selection reads candidate metadata and scores, never assay labels."""
 g=d.copy();g['_rank']=np.asarray(score);g=g.loc[mask] if mask is not None else g
 result=np.zeros(len(d),bool)
 for _,h in g.groupby('cluster',sort=True):
  n=int(np.floor(.25*len(h)));ix=h.sort_values(['_rank','row_id'],kind='stable').index[:n];result[ix]=True
 return result

def group_sum(d,columns,clusters):
 x=pd.DataFrame(columns,index=d.index);x['cluster']=d.cluster
 return x.groupby('cluster',sort=True).sum().reindex(clusters,fill_value=0).to_numpy(float)

class DatasetAudit:
 def __init__(self,d,name,draws=10000):
  self.d=d.reset_index(drop=True).copy();self.name=name;self.source=name=='Kircher_LOFO';self.formal=not self.source
  self.d['cluster']=self.d.fit_group if self.source else self.d.genomic_block
  self.clusters=sorted(self.d.cluster.unique());seed=int(hashlib.sha256(('913071'+name).encode()).hexdigest()[:8],16)
  self.w=weights(len(self.clusters),seed,draws);self.records=[];self.cache={};self.retention=[];self.point_metrics=[]
  self.err={(m,s):errors(self.d,m,s,self.source) for m in MODELS for s in ['all','strong','weak']}
  self.global_selection={}
  for m in MODELS:
   self.global_selection[m]=[select(self.d,self.d['risk_'+m]),select(self.d,-abs(self.d[m+'_DNASE_mean']))]
   for cl,h in self.d.groupby('cluster'):
    c,b=[int(mask[h.index].sum()) for mask in self.global_selection[m]]
    if c!=b or c!=int(np.floor(.25*len(h))):raise AssertionError('Global retention changed')
    self.retention.append({'dataset':name,'model':m,'mode':'global','cluster':cl,'candidate_count':len(h),'policy_retained':c,'magnitude_retained':b,'equal':c==b})
 def mask(self,model,dimension,level):
  if dimension=='all':return np.ones(len(self.d),bool)
  if dimension=='effect_strength':return np.isfinite(self.err[(model,level)])
  col=model+'_cell_match' if dimension=='cell_match' else dimension
  return self.d[col].eq(level).to_numpy()
 def sums(self,model,mode='global',dimension='all',level='all'):
  key=(model,mode,dimension,level)
  if key in self.cache:return self.cache[key]
  mask=self.mask(model,dimension,level)
  c,b=self.global_selection[model] if mode=='global' else [select(self.d,self.d['risk_'+model],mask),select(self.d,-abs(self.d[model+'_DNASE_mean']),mask)]
  err=self.err[(model,level if dimension=='effect_strength' else 'all')];valid=np.isfinite(err)&mask;ce=c&valid;be=b&valid
  x=group_sum(self.d,{'Cerror':np.where(ce,err,0),'Ccalls':ce.astype(int),'Merror':np.where(be,err,0),'Mcalls':be.astype(int),'Cretained':(c&mask).astype(int),'Mretained':(b&mask).astype(int),'population':mask.astype(int)},self.clusters)
  if mode=='metadata':
   if not np.array_equal(x[:,4],x[:,5]):raise AssertionError('Within-metadata retention differs')
   for i,cl in enumerate(self.clusters):self.retention.append({'dataset':self.name,'model':model,'mode':mode,'dimension':dimension,'level':level,'cluster':cl,'candidate_count':int(x[i,6]),'policy_retained':int(x[i,4]),'magnitude_retained':int(x[i,5]),'equal':x[i,4]==x[i,5]})
  self.cache[key]=x;return x
 def gain_record(self,m,mode,dim,level):
  x=self.sums(m,mode,dim,level);total=x.sum(axis=0);stat=float(gain(total));bs=gain(self.w@x[:,:4]);shared=int(((x[:,1]>0)&(x[:,3]>0)).sum());enough=total[1]>=20 and total[3]>=20 and shared>=10
  rec={'dataset':self.name,'model':m,'analysis_mode':mode,'dimension':dim,'level':level,'test_kind':'gain','test_id':'|'.join([self.name,m,mode,dim,level]),'formal_external':self.formal,**infer(stat,bs,enough),'population_variants':int(total[6]),'policy_retained':int(total[4]),'magnitude_retained':int(total[5]),'policy_audited_calls':int(total[1]),'magnitude_audited_calls':int(total[3]),'policy_sign_error':total[0]/total[1] if total[1] else np.nan,'magnitude_sign_error':total[2]/total[3] if total[3] else np.nan,'audit_opportunity_ratio':total[1]/total[3] if total[3] else np.nan,'contributing_clusters':shared,'orientation':'positive=confidence policy reduces error versus magnitude','status':'ESTIMABLE' if enough else 'NOT ESTIMABLE'}
  self.records.append(rec);return rec
 def interaction(self,m,mode,dim):
  a,b=[self.sums(m,mode,dim,lv) for lv in LEVELS[dim]];ta=a.sum(axis=0);tb=b.sum(axis=0);stat=float(gain(ta)-gain(tb));bs=gain(self.w@a[:,:4])-gain(self.w@b[:,:4])
  enough=all(x[:,1].sum()>=20 and x[:,3].sum()>=20 and ((x[:,1]>0)&(x[:,3]>0)).sum()>=10 for x in [a,b])
  self.records.append({'dataset':self.name,'model':m,'analysis_mode':mode,'dimension':dim,'level':'difference','test_kind':'interaction','test_id':'|'.join([self.name,m,mode,dim,'interaction']),'formal_external':self.formal,**infer(stat,bs,enough),'orientation':'gain('+LEVELS[dim][0]+') minus gain('+LEVELS[dim][1]+')','status':'ESTIMABLE' if enough else 'NOT ESTIMABLE'})
 def point_sums(self,m):
  x=self.d[m+'_DNASE_mean'].to_numpy(float);y=self.d.effect.to_numpy(float)
  return group_sum(self.d,{'n':np.ones(len(x)),'x':x,'y':y,'xx':x*x,'yy':y*y,'xy':x*y},self.clusters)
 def model_contrasts(self):
  pe,pb=[self.point_sums(m) for m in MODELS];qe=self.w@pe;qb=self.w@pb
  for m,ss,bs in [('enformer',pe,qe),('borzoi',pb,qb)]:
   values=pearson(bs);q=values[np.isfinite(values)];ci=np.quantile(q,[.025,.975]) if len(q) else [np.nan,np.nan]
   self.point_metrics.append({'dataset':self.name,'model':m,'readout':'DNASE','variants':len(self.d),'Pearson':float(pearson(ss.sum(axis=0))),'ci95_low':ci[0],'ci95_high':ci[1],'interval_type':'Pointwise estimation interval; no separate significance claim'})
  stat=float(pearson(pb.sum(axis=0))-pearson(pe.sum(axis=0)));bs=pearson(qb)-pearson(qe);enough=len(self.d)>=20 and len(self.clusters)>=10
  self.records.append({'dataset':self.name,'model':'paired','analysis_mode':'all_retention','dimension':'model','level':'Borzoi-Enformer Pearson','test_kind':'model_contrast','test_id':self.name+'|model|Pearson','formal_external':self.formal,**infer(stat,bs,enough),'orientation':'positive=Borzoi has higher Pearson','status':'ESTIMABLE' if enough else 'NOT ESTIMABLE'})
  e=self.err[('enformer','all')];b=self.err[('borzoi','all')];valid=np.isfinite(e)&np.isfinite(b)
  x=group_sum(self.d,{'errorE':np.where(valid,e,0),'nE':valid.astype(int),'errorB':np.where(valid,b,0),'nB':valid.astype(int)},self.clusters)
  total=x.sum(axis=0);st=total[0]/total[1]-total[2]/total[3] if total[1] else np.nan;draw=self.w@x;bs=np.divide(draw[:,0]-draw[:,2],draw[:,1],out=np.full(len(draw),np.nan),where=draw[:,1]>0)
  enough=valid.sum()>=20 and (x[:,1]>0).sum()>=10
  self.records.append({'dataset':self.name,'model':'paired','analysis_mode':'all_retention','dimension':'model','level':'Enformer-Borzoi sign error','test_kind':'model_contrast','test_id':self.name+'|model|sign','formal_external':self.formal,**infer(st,bs,enough),'paired_audited_variants':int(valid.sum()),'orientation':'positive=Borzoi has lower error on identical variants','status':'ESTIMABLE' if enough else 'NOT ESTIMABLE'})
  xe,xb=[self.sums(m) for m in MODELS];st=float(gain(xb.sum(axis=0))-gain(xe.sum(axis=0)));bs=gain(self.w@xb[:,:4])-gain(self.w@xe[:,:4]);enough=all(x[:,1].sum()>=20 and x[:,3].sum()>=20 and ((x[:,1]>0)&(x[:,3]>0)).sum()>=10 for x in [xe,xb])
  self.records.append({'dataset':self.name,'model':'paired','analysis_mode':'global','dimension':'model','level':'Borzoi-Enformer filtering gain','test_kind':'model_contrast','test_id':self.name+'|model|filtering','formal_external':self.formal,**infer(st,bs,enough),'orientation':'positive=Borzoi-specific policy has greater incremental gain','status':'ESTIMABLE' if enough else 'NOT ESTIMABLE'})
 def run(self):
  for m in MODELS:
   self.gain_record(m,'global','all','all')
   for dim in LEVELS:
    for lv in LEVELS[dim]:self.gain_record(m,'global',dim,lv)
    self.interaction(m,'global',dim)
   for dim in ['regulatory_type','cell_match']:
    for lv in LEVELS[dim]:self.gain_record(m,'metadata',dim,lv)
    self.interaction(m,'metadata',dim)
  self.model_contrasts();return self.records
