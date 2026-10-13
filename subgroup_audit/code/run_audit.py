"""CPU-only follow-up audit of completed Wahab scores. No new GPU inference."""
import argparse,json,os,hashlib,sys,time
from pathlib import Path
import numpy as np,pandas as pd,joblib
from scipy.stats import pearsonr,spearmanr
from common import ROOT,MODELS,MODS,plan,sha,write_json,input_reader
from annotate import annotate
from policies import fit,predict,source_lofo
from audit_engine import DatasetAudit
from audit_stats import holm,FAMILY
from noise_ceiling import audit as noise_audit,library_agreement

def verify_assets():
 p=ROOT/'config/analysis_assets_lock.json'
 if not p.exists():raise ValueError('Distribution lacks its analysis-asset lock')
 for name,digest in json.loads(p.read_text()).items():
  if sha(ROOT/name)!=digest:raise ValueError('Frozen analysis asset changed: '+name)

def modality_diagnostics(d):
 rows=[]
 for dataset,g in d.groupby('dataset'):
  conditions=list(g.groupby('condition')) if dataset=='Kircher' else [('all',g)]
  for cond,h in conditions:
   for model in MODELS:
    for mod in MODS:
     x=pd.to_numeric(h[model+'_'+mod+'_mean'],errors='coerce');mask=np.isfinite(x)&np.isfinite(pd.to_numeric(h.effect,errors='coerce'))
     valid=mask.sum()>=4 and x[mask].std()>0 and h.effect[mask].std()>0
     rows.append({'dataset':dataset,'condition':cond,'model':model,'modality':mod,'variants':int(mask.sum()),'Pearson':float(pearsonr(x[mask],h.effect[mask]).statistic) if valid else np.nan,'Spearman':float(spearmanr(x[mask],h.effect[mask]).statistic) if valid else np.nan,'primary_inferential_readout':mod=='DNASE','role':'Descriptive only; no oracle modality selection or uncorrected superiority claim'})
 return pd.DataFrame(rows)

def run(base,out=None):
 started=time.monotonic();p=plan();verify_assets();base=Path(base).resolve();out=Path(out).resolve() if out else ROOT/'results'
 if out==base or base in out.parents:raise ValueError('Outputs cannot be written inside the original score package')
 out.mkdir(parents=True,exist_ok=True)
 d,provenance=input_reader(base)
 identity={'plan_sha256':sha(ROOT/'config/preregistered_plan.json'),'assets_lock_sha256':sha(ROOT/'config/analysis_assets_lock.json'),'input_fingerprint':hashlib.sha256(json.dumps(provenance,sort_keys=True).encode()).hexdigest(),'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__}
 execution=out/'execution.json'
 if execution.exists() and json.loads(execution.read_text())['identity']!=identity:raise ValueError('Resume identity differs. Keep a separate result version.')
 if execution.exists() and json.loads(execution.read_text()).get('complete'):
  print('CPU audit already complete for this exact identity. Original FAIL unchanged.');return
 write_json({'identity':identity,'complete':False,'native_model_inference_performed':False,'input_provenance':provenance},execution)
 d,tracks=annotate(d,base);tracks.to_csv(out/'cell_match_provenance.csv',index=False)
 primary=d.loc[d.primary].copy();source=primary.loc[primary.dataset.eq('Kircher')]
 # New Borzoi policy is source-only. The old Enformer external policy is reused.
 bp=fit(source,'borzoi');joblib.dump(bp,out/'frozen_borzoi_source_policy.joblib')
 write_json({k:v for k,v in bp.items() if k not in ['estimator','calibrator']},out/'borzoi_source_policy_manifest.json')
 ext=primary.dataset.isin(['Abell','Ulirsch']);primary.loc[ext,'risk_enformer']=primary.loc[ext,'original_enformer_risk'];primary.loc[ext,'risk_borzoi']=predict(primary.loc[ext],bp)
 lofo={}
 for model in MODELS:
  risk,manifest=source_lofo(source,model);primary.loc[source.index,'risk_'+model]=risk;lofo[model]=manifest
 write_json(lofo,out/'source_LOFO_manifest.json')
 if not np.isfinite(primary[['risk_enformer','risk_borzoi']].to_numpy()).all():raise ValueError('Missing frozen policy ranks')
 allrecords=[];retained=[];points=[]
 for dataset in ['Abell','Ulirsch','Kircher']:
  name='Kircher_LOFO' if dataset=='Kircher' else dataset
  q=primary.loc[primary.dataset.eq(dataset)];a=DatasetAudit(q,name,p['bootstrap_draws']);a.run()
  allrecords.extend(a.records);retained.extend(a.retention);points.extend(a.point_metrics)
  print('Completed',name,'contrasts; source diagnostics cannot rescue the main FAIL.',flush=True)
  # A per-study result checkpoint survives an interrupted final write.
  write_json({'identity':identity,'dataset':name,'records':a.records,'point_metrics':a.point_metrics},out/('checkpoint_'+name+'.json'))
 z=pd.DataFrame(allrecords);formal=z.formal_external
 if int(formal.sum())!=FAMILY:raise AssertionError('Formal external family is not the preregistered70 tests')
 z['p_Holm']=np.nan;z.loc[formal,'p_Holm']=holm(z.loc[formal,'p_raw'].to_numpy());z['family_size']=np.where(formal,FAMILY,0)
 z['interpretation']='No resolved difference'
 z.loc[~z.evaluable,'interpretation']='NOT ESTIMABLE'
 z.loc[~formal,'interpretation']='Source LOFO diagnostic only'
 resolved=formal&z.evaluable&z.p_Holm.lt(.05)
 z.loc[resolved&z.test_kind.eq('gain')&z.estimate.lt(0),'interpretation']='Holm-resolved harm'
 benefit=resolved&z.test_kind.eq('gain')&z.estimate.gt(0)&z.audit_opportunity_ratio.ge(.9)
 z.loc[benefit,'interpretation']='Dataset-specific secondary benefit; independent replication required'
 z.loc[resolved&z.test_kind.ne('gain'),'interpretation']='Holm-resolved contrast; not a primary rescue'
 z.to_csv(out/'subgroup_tests.csv',index=False)
 pd.DataFrame(retained).to_csv(out/'retention_checks.csv',index=False);pd.DataFrame(points).to_csv(out/'point_model_estimates.csv',index=False)
 # Effects define diagnostic targets only; no validation value was fitted.
 primary.to_csv(out/'annotated_audit_variants.csv.gz',index=False,compression={'method':'gzip','mtime':0})
 noise=noise_audit(primary,p['bootstrap_draws']);noise.to_csv(out/'noise_ceiling_by_dataset.csv',index=False);library_agreement(d).to_csv(out/'Kircher_two_library_agreement.csv',index=False)
 modality_diagnostics(primary).to_csv(out/'descriptive_point_accuracy_by_modality.csv',index=False)
 coverage=[]
 for ds,g in primary.groupby('dataset'):
  for model in MODELS:
   for dim in ['regulatory_type','cell_match']:
    col=dim if dim=='regulatory_type' else model+'_cell_match'
    for level,n in g[col].value_counts(dropna=False).items():coverage.append({'dataset':ds,'model':model,'dimension':dim,'level':level,'variants':int(n)})
 pd.DataFrame(coverage).to_csv(out/'subgroup_population_coverage.csv',index=False)
 reconciliation=[]
 for dataset,old in provenance['original_primary_summary']['datasets'].items():
  r=z.loc[z.dataset.eq(dataset)&z.model.eq('enformer')&z.dimension.eq('all')&z.analysis_mode.eq('global')].iloc[0]
  if abs(r.estimate-old['primary_error_reduction'])>1e-10:raise AssertionError('Original Enformer primary point estimate changed')
  if int(r.policy_audited_calls)!=old['retained_audited_variants_candidate'] or int(r.magnitude_audited_calls)!=old['retained_audited_variants_baseline']:raise AssertionError('Original audited call counts changed')
  reconciliation.append({'dataset':dataset,'old_primary_gain':old['primary_error_reduction'],'recomputed_gain':r.estimate,'unchanged':True,'total_retained_policy':int(r.policy_retained),'total_retained_baseline':int(r.magnitude_retained),'audited_policy':int(r.policy_audited_calls),'audited_baseline':int(r.magnitude_audited_calls)})
 pd.DataFrame(reconciliation).to_csv(out/'original_primary_reconciliation.csv',index=False)
 from figures import build_figures
 build_figures(z,noise,pd.DataFrame(points),out)
 summary={'original_primary_decision':provenance['original_primary_summary']['decision'],'primary_verdict_changed':False,'follow_up_status':'COMPLETE','preregistered_external_family':FAMILY,'Holm_alpha':.05,'holm_resolved_gain_rows':int((resolved&z.test_kind.eq('gain')).sum()),'new_GPU_inference':False,'AlphaGenome_input_training_or_fitting':False,'noise_ceiling':noise.to_dict('records'),'source_inference':'Grouped LOFO diagnostics, conditional on overlapping source fits; not independent confirmatory validation','statement':'Subgroup discoveries are exploratory follow-up after a known global FAIL. No significant cell/type/effect subgroup can rescue the original primary result.'}
 write_json(summary,out/'follow_up_summary.json')
 (out/'SUMMARY.txt').write_text('Original primary: '+str(summary['original_primary_decision'])+' (unchanged). Subgroup follow-up COMPLETE.\nAll70 fixed external tests use Holm correction. Non-estimable strata remain in the correction family.\nTotal retained variants match; audited sign calls may differ.\nNo new GPU scoring or AlphaGenome outputs used.\n')
 write_json({f.name:sha(f) for f in sorted(out.glob('*.csv'))},out/'result_CSV_sha256.json')
 write_json({'identity':identity,'complete':True,'native_model_inference_performed':False,'wall_seconds':time.monotonic()-started,'input_provenance':provenance},execution)
 print('Original primary remains',summary['original_primary_decision'],'; CPU follow-up complete. No automatic positive paper claim.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--base',default=os.environ.get('RELIABILITY_BASE_DIR'));p.add_argument('--out');a=p.parse_args()
 if not a.base:p.error('Set RELIABILITY_BASE_DIR or pass --base /path/to/completed/Reliability_Audit_v2')
 run(a.base,a.out)
