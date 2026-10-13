"""Train on Kircher only, freeze, then audit two independent MPRA studies.

Deployment scores do not use assay spread, counts, effects or replicate labels.
Replicates define audit targets only. Confidence intervals resample 1 Mb blocks,
not nearby variants or multiple RNA views as independent biological samples.
"""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr,spearmanr,rankdata
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.isotonic import IsotonicRegression
import joblib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import ROOT,SEED,load_protocol,sha,atomic_json,permission_for_alpha,code_hashes
METHODS=['Source-trained risk','Magnitude only','Cross-track disagreement','Cross-model disagreement','Random ranking']
def clean_json(v):
    if isinstance(v,dict):return {k:clean_json(x) for k,x in v.items()}
    if isinstance(v,list):return [clean_json(x) for x in v]
    if isinstance(v,np.generic):v=v.item()
    if isinstance(v,float) and not np.isfinite(v):return None
    return v
def collect(root=ROOT,shards=8):
    d=pd.read_csv(root/'data/variants.csv');frames=[]
    for model in ['enformer','borzoi']:
        rows=[]
        for s in range(shards):
            state=root/'scores'/model/f'shard{s:02d}.json'
            if not state.exists() or not json.loads(state.read_text()).get('complete'):raise RuntimeError('INCOMPLETE '+model+' shard '+str(s)+'. Resume before evaluating.')
            status=json.loads(state.read_text())
            if status['identity']['source_rows_sha256']!=sha(root/'data/variants.csv'):raise ValueError('Scored row universe changed')
            if status['identity']['protocol_sha256']!=sha(root/'config/protocol.json'):raise ValueError('Scored protocol changed')
            if status['identity']['code_sha256']!=code_hashes():raise ValueError('Scoring implementation changed')
            if not status.get('real_model_scores'):raise ValueError('Non-model/test outputs prohibited in scientific audit')
            rows.extend(json.loads(x) for x in (root/'scores'/model/f'shard{s:02d}.jsonl').read_text().splitlines())
        z=pd.DataFrame(rows)
        if not z.row_id.is_unique or set(z.row_id)!=set(d.row_id):raise ValueError('Missing/duplicate model scores')
        if not z.real_model_scores.eq(True).all():raise ValueError('Synthetic outputs cannot produce a scientific verdict')
        z=z.drop(columns=['model','real_model_scores','variant_id','scoring_device','cuda_version']).rename(columns={c:model+'_'+c for c in z if c not in ['row_id','variant_id','model','real_model_scores','scoring_device','cuda_version']})
        frames.append(z)
    return d.merge(frames[0],on='row_id',validate='one_to_one').merge(frames[1],on='row_id',validate='one_to_one')
def deployment_features(d,scale=None,alpha=False):
    """Every input is a model output; no observed assay quantity enters."""
    e=d.enformer_DNASE_mean.to_numpy(float);b=d.borzoi_DNASE_mean.to_numpy(float)
    if not np.isfinite(e).all() or not np.isfinite(b).all():raise ValueError('Primary DNASE scores missing')
    if scale is None:
        scale={'e':max(1e-6,float(np.std(e))),'b':max(1e-6,float(np.std(b)))}
    track=d.enformer_DNASE_std.to_numpy(float)/scale['e']
    cross=abs(e/scale['e']-b/scale['b'])
    # Context spread is label-free but candidate-set-dependent. Its unit is a
    # fixed genomic block, never a grouping selected using observed effects.
    keys=[d.dataset,d.genomic_block]
    espal=d.enformer_DNASE_mean.groupby(keys).transform(lambda v:float(np.std(v,ddof=0))).to_numpy()/scale['e']
    bspread=d.borzoi_DNASE_mean.groupby(keys).transform(lambda v:float(np.std(v,ddof=0))).to_numpy()/scale['b']
    # Fixed signed modality summaries capture agreement, with explicit missing
    # indicators instead of pretending an absent assay is a zero prediction.
    cols=[]
    for model in ['enformer','borzoi']:
        for mod in ['ATAC','CAGE','H3K27ac']:
            v=d[model+'_'+mod+'_mean'].to_numpy(float)
            cols.extend([np.nan_to_num(v,nan=0),np.isfinite(v).astype(float)])
    if alpha:
        permission_for_alpha()
        ag=d.alphagenome_DNASE_mean.to_numpy(float)
        if not np.isfinite(ag).all():raise ValueError('Incomplete optional AlphaGenome population')
        if 'ag' not in scale:scale['ag']=max(1e-6,float(np.std(ag)))
        ac=ag/scale['ag'];cross=np.maximum.reduce([cross,abs(e/scale['e']-ac),abs(b/scale['b']-ac)])
        cols.extend([ac])
    x=np.column_stack([np.log1p(abs(e)/scale['e']),track,cross,espal,bspread,*cols])
    return x,scale,{'magnitude':abs(e),'track':track,'cross':cross,'forecast':e}
def balanced(groups):
    counts=pd.Series(groups).value_counts();w=np.array([1/counts[g] for g in groups]);return w/w.mean()
def fit_policy(d,root=ROOT,alpha=False):
    if alpha:permission_for_alpha()
    p=load_protocol(root);source=d.loc[d.dataset.eq('Kircher') & d.primary].copy()
    groups=sorted(source.fit_group.unique(),key=lambda s:hashlib.sha256((str(SEED)+s).encode()).hexdigest())
    calgroups=groups[:int(np.ceil(.25*len(groups)))];fitgroups=groups[len(calgroups):]
    fit=source.loc[source.fit_group.isin(fitgroups)];cal=source.loc[source.fit_group.isin(calgroups)]
    # Prediction magnitude scales use fit-family outputs only, with no labels.
    _,scale,_=deployment_features(fit,alpha=alpha);x,_,aux=deployment_features(fit,scale,alpha=alpha)
    q=fit.effect.abs().ge(.5).to_numpy();y=(np.sign(aux['forecast'][q])!=np.sign(fit.effect.to_numpy()[q])).astype(int)
    if q.sum()<20:raise ValueError('Insufficient Kircher sign training observations')
    model=HistGradientBoostingClassifier(max_iter=120,max_leaf_nodes=7,min_samples_leaf=30,l2_regularization=5,learning_rate=.07,random_state=SEED)
    constant=float(y.mean()) if len(np.unique(y))<2 else None
    if constant is None:model.fit(x[q],y,sample_weight=balanced(fit.fit_group.to_numpy()[q]))
    z,_,a=deployment_features(cal,scale,alpha=alpha);mask=cal.effect.abs().ge(.5).to_numpy()
    raw=np.full(len(z),constant) if constant is not None else model.predict_proba(z)[:,1]
    cy=(np.sign(a['forecast'][mask])!=np.sign(cal.effect.to_numpy()[mask])).astype(int)
    iso=None
    if mask.sum() and np.unique(raw[mask]).size>1:
        iso=IsotonicRegression(out_of_bounds='clip').fit(raw[mask],cy,sample_weight=balanced(cal.fit_group.to_numpy()[mask]))
    # Source-frozen bins are used only by the conditional magnitude audit.
    boundaries=np.unique(np.quantile(abs(fit.enformer_DNASE_mean),[.2,.4,.6,.8])).tolist()
    policy={'model':model,'constant':constant,'iso':iso,'scale':scale,'magnitude_boundaries':boundaries,'fit_groups':fitgroups,'calibration_groups':calgroups,'fit_rows':len(fit),'fit_effectful_rows':int(q.sum()),'calibration_effectful_rows':int(mask.sum()),'protocol_sha256':sha(root/'config/protocol.json'),'source_prediction_hash':hashlib.sha256(source.drop(columns=[c for c in source if c.startswith('audit_')]).to_csv(index=False).encode()).hexdigest(),'uses_AlphaGenome':alpha,'feature_kind':'Model-output summaries only; measured assay spread, counts and replicate information prohibited as features'}
    (root/'models').mkdir(exist_ok=True);tag='alpha_optional' if alpha else 'open'
    joblib.dump(policy,root/('models/frozen_'+tag+'_policy.joblib'))
    atomic_json({k:v for k,v in policy.items() if k not in ['model','iso']},root/('models/frozen_'+tag+'_policy_manifest.json'))
    return policy
def apply_policy(d,policy):
    x,_,a=deployment_features(d,policy['scale'],alpha=policy.get('uses_AlphaGenome',False))
    raw=np.full(len(x),policy['constant']) if policy['constant'] is not None else policy['model'].predict_proba(x)[:,1]
    prob=policy['iso'].predict(raw) if policy['iso'] is not None else raw
    result=d.copy();result['forecast']=a['forecast'];result['risk_raw']=raw;result['risk_probability']=prob
    result['score|Source-trained risk']=raw;result['score|Magnitude only']=-a['magnitude'];result['score|Cross-track disagreement']=a['track'];result['score|Cross-model disagreement']=a['cross']
    result['score|Random ranking']=result.row_id.map(lambda s:int(hashlib.sha256((str(SEED)+s).encode()).hexdigest()[:12],16)/16**12)
    result['magnitude_bin']=np.searchsorted(policy['magnitude_boundaries'],a['magnitude'],side='right')
    # Screening and sign truth use disjoint replicate samples. Each variant is
    # counted once, with its error averaged over its eligible audit views.
    errors=[];views=[]
    for r in result.itertuples():
        vals=[]
        for label in 'ABC':
            screen=getattr(r,'audit_screen_'+label);truth=getattr(r,'audit_eval_'+label)
            if np.isfinite(screen) and np.isfinite(truth) and abs(screen)>=.5 and truth!=0:vals.append(float(np.sign(r.forecast)!=np.sign(truth)))
        errors.append(float(np.mean(vals)) if vals else np.nan);views.append(len(vals))
    result['audit_sign_error']=errors;result['audit_views']=views
    return result
def select(g,method,rate):
    n=int(np.floor(rate*len(g)))
    return g.sort_values(['score|'+method,'row_id'],kind='stable').index[:n]
def partial_rank(x,y,controls):
    z=np.column_stack([np.ones(len(x)),*[rankdata(c) for c in controls]])
    rx=rankdata(x);ry=rankdata(y)
    a=rx-z@np.linalg.lstsq(z,rx,rcond=None)[0];b=ry-z@np.linalg.lstsq(z,ry,rcond=None)[0]
    if len(a)<5 or np.std(a)<1e-12 or np.std(b)<1e-12:return np.nan
    return float(pearsonr(a,b).statistic)
def noise_corrected(g):
    if len(g)<10:return np.nan,np.nan,np.nan,'fewer than10 variants'
    y=g.effect.to_numpy(float);p=g.forecast.to_numpy(float)
    var=float(np.var(y,ddof=1));noise=float(g.noise_variance.mean())
    rel=1-noise/var if var>0 else np.nan
    raw=float(pearsonr(y,p).statistic) if np.std(p)>0 and np.std(y)>0 else np.nan
    if not np.isfinite(rel) or rel<=0:return raw,rel,np.nan,'noise ceiling unidentified'
    corrected=raw/np.sqrt(rel)
    return raw,rel,corrected,'unstable corrected correlation >1' if abs(corrected)>1 else 'classical additive-error approximation'
def summarize_selection(g,idx):
    p=g.loc[idx];a=p.loc[p.audit_sign_error.notna()]
    return {'total':len(g),'retained':len(p),'audit_calls':len(a),'error_sum':float(a.audit_sign_error.sum()),'sign_error':float(a.audit_sign_error.mean()) if len(a) else np.nan,'sign_accuracy':1-float(a.audit_sign_error.mean()) if len(a) else np.nan,'mean_held_out_views':float(a.audit_views.mean()) if len(a) else np.nan}
def paired_bootstrap(candidate,baseline,rng,n=10000):
    """Aggregate sign risk, with paired block resampling and variable call Ns."""
    q=candidate.merge(baseline,on='block',suffixes=('_C','_B'),validate='one_to_one')
    z=q[['error_sum_C','audit_calls_C','error_sum_B','audit_calls_B']].to_numpy(float)
    draws=rng.integers(0,len(z),size=(n,len(z)));s=z[draws].sum(axis=1)
    good=(s[:,1]>0)&(s[:,3]>0);delta=s[good,2]/s[good,3]-s[good,0]/s[good,1]
    point=z[:,2].sum()/z[:,3].sum()-z[:,0].sum()/z[:,1].sum() if z[:,1].sum()>0 and z[:,3].sum()>0 else np.nan
    return point,np.quantile(delta,[.025,.975]).tolist() if len(delta) else [np.nan,np.nan],np.quantile(delta,[.0125,.9875]).tolist() if len(delta) else [np.nan,np.nan]
def evaluate(d,policy,root=ROOT,boot=10000,output_dir=None):
    p=load_protocol(root);out=Path(output_dir) if output_dir else root/'results';out.mkdir(parents=True,exist_ok=True);rng=np.random.default_rng(812071)
    ext=d.loc[d.dataset.isin(['Ulirsch','Abell']) & d.primary].copy()
    pred=apply_policy(ext,policy);pred.to_csv(out/'external_predictions.csv',index=False)
    metrics=[];conditional=[];units=[];comparisons=[];gates={}
    for (dataset,block),g in pred.groupby(['dataset','genomic_block'],sort=True):
        for rate in [.10,.25,.50,.75]:
            for method in METHODS:
                record={'dataset':dataset,'block':block,'retention':rate,'method':method,**summarize_selection(g,select(g,method,rate))};metrics.append(record)
        # Same perblock and permagnitude-bin budget as the magnitude baseline.
        # The budget depends on predictions only, not labels or measured effects.
        base=select(g,'Magnitude only',.25);chosen=[]
        for b,h in g.groupby('magnitude_bin'):
            n=int(g.loc[base,'magnitude_bin'].eq(b).sum())
            chosen.extend(h.sort_values(['score|Source-trained risk','row_id']).index[:n])
        conditional.append({'dataset':dataset,'block':block,'retention':.25,'method':'Conditional risk',**summarize_selection(g,chosen)})
        raw,rel,corr,status=noise_corrected(g)
        units.append({'dataset':dataset,'block':block,'variants':len(g),'observed_correlation':raw,'mean_estimate_reliability':rel,'noise_corrected_correlation':corr,'noise_status':status,'measured_spread':g.effect.std(),'mean_noise_SE':np.sqrt(g.noise_variance.mean()),'prediction_spread':g.forecast.std(),'borzoi_prediction_spread':g.borzoi_DNASE_mean.std(),'cross_track_disagreement':g['score|Cross-track disagreement'].mean(),'cross_model_disagreement':g['score|Cross-model disagreement'].mean()})
    m=pd.DataFrame(metrics);cond=pd.DataFrame(conditional);unit=pd.DataFrame(units)
    m.to_csv(out/'matched_retention_by_block.csv',index=False);cond.to_csv(out/'conditional_magnitude_selection.csv',index=False);unit.to_csv(out/'variance_robust_accuracy_by_block.csv',index=False)
    if len(m) and not m.groupby(['dataset','block','retention']).retained.nunique().eq(1).all():raise AssertionError('Unequal retention')
    for dataset,g in pred.groupby('dataset',sort=True):
        z=m.loc[m.dataset.eq(dataset)&m.retention.eq(.25)];base=z.loc[z.method.eq('Magnitude only')];risk=z.loc[z.method.eq('Source-trained risk')]
        delta,ci,sim=paired_bootstrap(risk,base,rng,boot)
        cd,cci,_=paired_bootstrap(cond.loc[cond.dataset.eq(dataset)],base,rng,boot)
        paired=risk.merge(base,on='block',suffixes=('_C','_B'))
        validblocks=int(((paired.audit_calls_C>0)&(paired.audit_calls_B>0)).sum())
        opportunities=risk.audit_calls.sum()/base.audit_calls.sum() if base.audit_calls.sum()>0 else 0
        # Allobserved-point fidelity is a diagnostic independent of retention.
        rows=[]
        # Block sufficient statistics give exactly the pooled weighted Pearson
        # correlation without O(N) mapping in every bootstrap iteration.
        ss=[]
        for _,h in g.groupby('genomic_block'):
            x=h.forecast.to_numpy(float);y=h.effect.to_numpy(float)
            ss.append([len(h),x.sum(),y.sum(),(x*x).sum(),(y*y).sum(),(x*y).sum()])
        ss=np.asarray(ss,float)
        for _ in range(boot):
            n,sx,sy,sxx,syy,sxy=ss[rng.integers(0,len(ss),len(ss))].sum(axis=0)
            denom=np.sqrt(max(0,sxx-sx*sx/n)*max(0,syy-sy*sy/n))
            if denom>0:rows.append((sxy-sx*sy/n)/denom)
        pclo,pchi=np.quantile(rows,[.025,.975]) if rows else [np.nan,np.nan]
        enough=(risk.audit_calls.sum()>=20 and base.audit_calls.sum()>=20 and validblocks>=10)
        success=bool(enough and sim[0]>0 and cci[0]>0 and opportunities>=.9 and pclo>0)
        gates[dataset]={'evaluability':bool(enough),'retained_audited_variants_candidate':int(risk.audit_calls.sum()),'retained_audited_variants_baseline':int(base.audit_calls.sum()),'contributing_blocks':validblocks,'primary_error_reduction':delta,'CI95':ci,'CI97_5':sim,'conditional_error_reduction':cd,'conditional_CI95':cci,'audit_opportunity_ratio':opportunities,'point_correlation_CI95':[pclo,pchi],'pass':success}
        comparisons.append({'dataset':dataset,'error_reduction':delta,'ci95_low':ci[0],'ci95_high':ci[1],'ci97_5_low':sim[0],'ci97_5_high':sim[1],'conditional_gain':cd,'conditional_ci95_low':cci[0],'conditional_ci95_high':cci[1],'retention':.25})
    pd.DataFrame(comparisons).to_csv(out/'primary_bootstrap_comparisons.csv',index=False)
    # Absolute precision and effect fidelity are reported, not just differences.
    absolute=[];points=[];bins=[]
    for dataset,g in pred.groupby('dataset'):
        for model in ['enformer','borzoi']:
            for mod in ['DNASE','ATAC','CAGE','H3K27ac']:
                h=g.loc[g[model+'_'+mod+'_mean'].notna()].copy();h['forecast']=h[model+'_'+mod+'_mean']
                raw,rel,corr,nstatus=noise_corrected(h)
                points.append({'dataset':dataset,'model':model,'modality':mod,'variants':len(h),'Pearson':raw,'Spearman':float(spearmanr(h.forecast,h.effect).statistic) if len(h)>3 and h.forecast.std()>0 else np.nan,'estimated_mean_reliability':rel,'noise_ceiling_corrected_Pearson':corr,'noise_status':nstatus,'primary_point_readout':model=='enformer' and mod=='DNASE'})
        for rate in [.10,.25,.50,.75]:
            for method in METHODS:
                q=m.loc[m.dataset.eq(dataset)&m.retention.eq(rate)&m.method.eq(method)]
                arr=q[['error_sum','audit_calls']].to_numpy(float);draw=arr[rng.integers(0,len(arr),(boot,len(arr)))].sum(axis=1);ok=draw[:,1]>0
                prec=1-draw[ok,0]/draw[ok,1];lo,hi=np.quantile(prec,[.025,.975]) if len(prec) else [np.nan,np.nan]
                total=arr.sum(axis=0);absolute.append({'dataset':dataset,'method':method,'retention':rate,'retained':int(q.retained.sum()),'audited_variants':int(total[1]),'precision':1-total[0]/total[1] if total[1] else np.nan,'ci95_low':lo,'ci95_high':hi})
        b=g.loc[g.audit_sign_error.notna()].copy();b['probability_bin']=pd.cut(b.risk_probability,bins=np.linspace(0,1,6),include_lowest=True)
        for band,h in b.groupby('probability_bin',observed=True):
            bins.append({'dataset':dataset,'bin':str(band),'variants':len(h),'predicted_error_probability':h.risk_probability.mean(),'observed_replicate_separated_error':h.audit_sign_error.mean()})
    pd.DataFrame(absolute).to_csv(out/'precision_with_bootstrap_CIs.csv',index=False)
    pd.DataFrame(points).to_csv(out/'point_accuracy_by_dataset.csv',index=False)
    pd.DataFrame(bins).to_csv(out/'probability_calibration.csv',index=False)
    adjusted=[]
    for dataset,g in unit.groupby('dataset',sort=True):
        q=g.replace([np.inf,-np.inf],np.nan).dropna(subset=['noise_corrected_correlation','measured_spread','mean_noise_SE'])
        q=q.loc[q.noise_corrected_correlation.abs().le(1)]
        for score in ['prediction_spread','borzoi_prediction_spread','cross_track_disagreement','cross_model_disagreement']:
            a=q.dropna(subset=[score]);rho=partial_rank(a[score],a.noise_corrected_correlation,[a.measured_spread,a.mean_noise_SE]) if len(a)>=8 else np.nan
            bs=[]
            if len(a)>=8:
                for _ in range(boot):
                    ix=rng.integers(0,len(a),len(a));h=a.iloc[ix];v=partial_rank(h[score],h.noise_corrected_correlation,[h.measured_spread,h.mean_noise_SE])
                    if np.isfinite(v):bs.append(v)
            lo,hi=np.quantile(bs,[.025,.975]) if bs else [np.nan,np.nan]
            adjusted.append({'dataset':dataset,'score':score,'audit_units':len(a),'partial_rank_correlation':rho,'ci95_low':lo,'ci95_high':hi,'controls':'Measured spread and replicate-derived mean SE; audit only, never inference features','note':'Exploratory block-level conditional association, not a learned deployment model'})
    pd.DataFrame(adjusted).to_csv(out/'partial_correlations.csv',index=False)
    complete=set(gates)==set(['Ulirsch','Abell'])
    decision='INCONCLUSIVE' if not complete or not all(x['evaluability'] for x in gates.values()) else 'PASS' if all(x['pass'] for x in gates.values()) else 'FAIL'
    status=clean_json({'decision':decision,'primary_model':'Enformer DNASE frozen point forecast','primary_policy':'Kircher-trained open-model risk','datasets':gates,'real_model_scores':True,'protocol_sha256':sha(root/'config/protocol.json'),'AlphaGenome_in_primary_fit':policy.get('uses_AlphaGenome',False),'analysis_role':'Exploratory permission-gated secondary; cannot rescue primary verdict' if policy.get('uses_AlphaGenome',False) else 'Primary locked open-model audit','clinical_interpretation':False,'bootstrap':'Paired 1 Mb block resampling, conditional on fixed policies; genomic blocks are not proven LD-independent','prospective_status':'Frozen on Kircher; no validation outcomes fitted to deployment policy','noise_correction_status':'Classical additive-error approximation; no pervariant SEs manufactured for Kircher'})
    atomic_json(status,out/'audit_summary.json')
    atomic_json({f.name:sha(f) for f in sorted(out.glob('*.csv'))},out/'result_CSV_sha256.json')
    (out/'SUMMARY.txt').write_text(decision+': '+('External improvement beyond magnitude passed in both datasets.' if decision=='PASS' else 'No reliable incremental benefit under the locked rule.' if decision=='FAIL' else 'Insufficient complete/evaluable evidence; do not convert this into a success or failure.')+'\n\n'+json.dumps(status,indent=2)+'\n')
    figures(pd.DataFrame(comparisons),pd.DataFrame(adjusted),m,out)
    if bins:
        f,a=plt.subplots(figsize=(4.8,4));a.plot([0,1],[0,1],'--',color='gray')
        bc=pd.DataFrame(bins)
        for ds,h in bc.groupby('dataset'):a.plot(h.predicted_error_probability,h.observed_replicate_separated_error,'o-',label=ds)
        a.set(xlabel='Source-calibrated error probability',ylabel='Observed sign-error fraction',xlim=(0,1),ylim=(0,1),title='External calibration diagnostic');a.legend();f.tight_layout()
        for suffix in ['png','pdf','svg']:f.savefig(out/('figure4_probability_calibration.'+suffix),dpi=320,bbox_inches='tight')
        plt.close(f)
    return status
def figures(corr,adjusted,metrics,out):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.right':False,'axes.spines.top':False,'savefig.dpi':320})
    def save(f,n):
        f.tight_layout()
        for e in ['png','pdf','svg']:f.savefig(out/(n+'.'+e),bbox_inches='tight')
        plt.close(f)
    f,a=plt.subplots(figsize=(6,3.5));v=corr.error_reduction.to_numpy()*100;lo=corr.ci97_5_low.to_numpy()*100;hi=corr.ci97_5_high.to_numpy()*100
    a.errorbar(v,np.arange(len(v)),xerr=np.vstack([np.maximum(0,v-lo),np.maximum(0,hi-v)]),fmt='o',capsize=4,color='#076A73');a.set_yticks(np.arange(len(v)),corr.dataset);a.axvline(0,ls='--',color='gray');a.set_xlabel('Sign-error reduction versus magnitude (percentage points)');a.set_title('25% retention; 97.5% paired block intervals');save(f,'figure1_external_gain')
    f,a=plt.subplots(figsize=(6.5,3.5))
    for dataset,g in metrics.groupby('dataset'):
        for method in ['Source-trained risk','Magnitude only']:
            q=g.loc[g.method.eq(method)].groupby('retention')[['error_sum','audit_calls']].sum();a.plot(q.index,1-q.error_sum/q.audit_calls,'o-',label=dataset+' / '+method)
    a.set(xlabel='Fraction retained within each genomic block',ylabel='Replicate-separated sign accuracy',title='Fixed effect-screen threshold 0.5 log2');a.legend(fontsize=7);save(f,'figure2_risk_retention')
    f,a=plt.subplots(figsize=(7,4));z=adjusted.replace([np.inf,-np.inf],np.nan).dropna(subset=['partial_rank_correlation'])
    v=z.partial_rank_correlation.to_numpy();lo=z.ci95_low.to_numpy();hi=z.ci95_high.to_numpy()
    a.errorbar(v,np.arange(len(z)),xerr=np.vstack([np.maximum(0,v-lo),np.maximum(0,hi-v)]),fmt='o',capsize=3);a.set_yticks(np.arange(len(z)),z.dataset+' / '+z.score.str.replace('_',' '));a.axvline(0,color='gray',ls='--');a.set_xlabel('Partial rank correlation after spread/noise adjustment');save(f,'figure3_adjusted_reliability')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--shards',type=int,default=8);p.add_argument('--bootstrap',type=int,default=10000);a=p.parse_args();d=collect(shards=a.shards)
    sources=[]
    for condition,g in d.loc[d.dataset.eq('Kircher') & d.primary].groupby('condition'):
        for model in ['enformer','borzoi']:
            for mod in ['DNASE','ATAC','CAGE','H3K27ac']:
                x=g[model+'_'+mod+'_mean'];mask=x.notna()
                sources.append({'condition':condition,'model':model,'modality':mod,'variants':int(mask.sum()),'Pearson':pearsonr(x[mask],g.effect[mask]).statistic if mask.sum()>3 and x[mask].std()>0 else np.nan,'Spearman':spearmanr(x[mask],g.effect[mask]).statistic if mask.sum()>3 and x[mask].std()>0 else np.nan,'primary_modality':mod=='DNASE','note':'Kircher source diagnostic only; not independent validation'})
    pd.DataFrame(sources).to_csv(ROOT/'results/source_point_accuracy.csv',index=False)
    policy=fit_policy(d);evaluate(d,policy,boot=a.bootstrap)
