"""Paired genomic-cluster bootstrap and fixed-family multiplicity correction.

Bootstrap tests are null-centered, two-sided approximations conditional on
fixed policies. Percentile95% intervals are pointwise; Bonferroni intervals
use the entire70-test external family. Neither is an exact dependent-data guarantee.
"""
import numpy as np
FAMILY=70

def weights(n,seed,draws=10000):
 if n<1:raise ValueError('No audit clusters')
 return np.random.default_rng(seed).multinomial(n,np.ones(n)/n,size=draws).astype(np.float64)

def ratio(e,n):
 return np.divide(e,n,out=np.full_like(np.asarray(e,dtype=float),np.nan),where=np.asarray(n)>0)

def gain(a):
 # Array columns: candidate errors, candidate calls, magnitude errors/calls.
 return ratio(a[...,2],a[...,3])-ratio(a[...,0],a[...,1])

def pearson(s):
 # n, sumX, sumY, sumX², sumY², sumXY
 n=s[...,0];numer=s[...,5]-s[...,1]*s[...,2]/n
 den=np.sqrt(np.maximum(0,s[...,3]-s[...,1]**2/n)*np.maximum(0,s[...,4]-s[...,2]**2/n))
 return ratio(numer,den)

def infer(point,draws,evaluable=True):
 v=np.asarray(draws,float);v=v[np.isfinite(v)];valid=len(v)/len(draws) if len(draws) else 0
 enough=bool(evaluable and np.isfinite(point) and valid>=.95)
 if not len(v):return {'estimate':point,'ci95_low':np.nan,'ci95_high':np.nan,'ci_family_low':np.nan,'ci_family_high':np.nan,'p_raw':1.0,'evaluable':False,'bootstrap_valid_fraction':valid}
 ci=np.quantile(v,[.025,.975]);bc=np.quantile(v,[.05/(2*FAMILY),1-.05/(2*FAMILY)])
 # Recenter the paired cluster-bootstrap distribution under the null. This is
 # an asymptotic test, not a permutation test or an exact randomization claim.
 p=(1+int((abs(v-point)>=abs(point)).sum()))/(len(v)+1) if enough else 1.0
 return {'estimate':point,'ci95_low':ci[0],'ci95_high':ci[1],'ci_family_low':bc[0],'ci_family_high':bc[1],'p_raw':p,'evaluable':enough,'bootstrap_valid_fraction':valid}

def holm(p):
 a=np.asarray(p,float);a=np.where(np.isfinite(a),a,1);n=len(a);order=np.argsort(a,kind='stable');out=np.empty(n);running=0
 for k,i in enumerate(order):running=max(running,(n-k)*a[i]);out[i]=min(1,running)
 return out

def noise(s):
 # n, sumY, sumY², sumSEM². Invalid noise data are handled before this call.
 n=s[...,0];var=ratio(s[...,2]-s[...,1]**2/n,n-1);mean_noise=ratio(s[...,3],n)
 rel=1-ratio(mean_noise,var);ceiling=np.sqrt(np.where(rel>0,rel,np.nan))
 return rel,ceiling
