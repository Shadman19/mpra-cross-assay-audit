import numpy as np, pandas as pd
R='/home/spath004/cardiac_r21/repo/mpra-cross-assay-audit'
d=pd.read_csv(R+'/subgroup_audit/results/annotated_audit_variants.csv.gz',low_memory=False)
d=d[d.primary.astype(str).str.lower().eq('true')].copy()
for c in ['effect','enformer_DNASE_mean','borzoi_DNASE_mean']:
    if c in d: d[c]=pd.to_numeric(d[c],errors='coerce')
v=pd.read_csv('/home/spath004/cardiac_r21/probe_results/malinois/malinois_variant_scores.csv.gz')
d=d.merge(v[['variant_id','malinois_K562','malinois_HepG2','malinois_SKNSH']],on='variant_id',how='left')
d['cluster']=np.where(d.dataset.eq('Kircher'),d.fit_group.astype(str),d.genomic_block.astype(str))
def sums(p,y,inv,nu):
    s=np.zeros((nu,6)); np.add.at(s,inv,np.column_stack([np.ones_like(y),p,y,p*p,y*y,p*y])); return s
def r_of(T):
    n=T[:,0]; return (T[:,5]-T[:,1]*T[:,2]/n)/np.sqrt((T[:,3]-T[:,1]**2/n)*(T[:,4]-T[:,2]**2/n))
for ds in ['Abell','Ulirsch']:
    t=d[d.dataset.eq(ds)]
    for col in ['malinois_K562','malinois_HepG2','malinois_SKNSH']:
        m=np.isfinite(t[col])&np.isfinite(t.enformer_DNASE_mean)&np.isfinite(t.effect)
        g=t[m]; y=g.effect.to_numpy(); u,inv=np.unique(g.cluster.to_numpy(),return_inverse=True)
        A=sums(g[col].to_numpy(),y,inv,len(u)); B=sums(g.enformer_DNASE_mean.to_numpy(),y,inv,len(u))
        W=np.random.default_rng(20261015).multinomial(len(u),np.ones(len(u))/len(u),size=4000)
        dif=r_of(W@A)-r_of(W@B); lo,hi=np.nanpercentile(dif,[2.5,97.5])
        pt=np.corrcoef(g[col],y)[0,1]-np.corrcoef(g.enformer_DNASE_mean,y)[0,1]
        print(f'{ds:8s} {col:16s} diff {pt:+.3f} [{lo:+.3f}, {hi:+.3f}]  n={len(g)} blocks={len(u)}')
