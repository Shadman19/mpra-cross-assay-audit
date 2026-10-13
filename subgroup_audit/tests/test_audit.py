"""Synthetic software checks. No fake outcomes are written into results/."""
import unittest,tempfile,sys,json
from pathlib import Path
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'code'))
from common import plan,sha
from audit_stats import holm,infer,noise,weights
from audit_engine import DatasetAudit,errors,select
from policies import fit,predict,source_lofo,features
from annotate import expected_cell

class AuditTests(unittest.TestCase):
 @staticmethod
 def data(dataset='Abell',groups=12,n=24):
  rng=np.random.default_rng(123);rows=[]
  for g in range(groups):
   for j in range(n):
    y=rng.normal(0,1);e=y+rng.normal(0,.5);b=y+rng.normal(0,.4)
    r={'row_id':dataset+'|'+str(g)+'|'+str(j),'dataset':dataset,'primary':True,'effect':y,'noise_variance':.05,'fit_group':'family'+str(g),'genomic_block':'chr1:'+str(g),'regulatory_type':'promoter' if j%2 else 'enhancer','enformer_cell_match':'matched','borzoi_cell_match':'matched','risk_enformer':abs(rng.normal()),'risk_borzoi':abs(rng.normal()),'enformer_DNASE_mean':e,'borzoi_DNASE_mean':b,'enformer_DNASE_std':.1,'borzoi_DNASE_std':.1}
    for m,p in [('enformer',e),('borzoi',b)]:
     for mod in ['ATAC','CAGE','H3K27ac']:r[m+'_'+mod+'_mean']=p+rng.normal(0,.1)
    for lab in 'ABC':r['audit_screen_'+lab]=y+rng.normal(0,.1);r['audit_eval_'+lab]=y+rng.normal(0,.1)
    rows.append(r)
  return pd.DataFrame(rows)
 def test_plan_locked_and_family_70(self):
  p=plan();self.assertEqual(p['formal_test_family']['total_fixed_tests'],70)
 def test_holm_known_values_includes_unestimable(self):
  np.testing.assert_allclose(holm([.001,.02,.3,1]),[.004,.06,.6,1])
 def test_null_identity_bootstrap(self):
  z=infer(0,np.zeros(10000),True);self.assertEqual(z['p_raw'],1);self.assertEqual(z['ci95_low'],0)
 def test_ceiling_no_clipping(self):
  r,c=noise(np.array([10,0,9,0.]));self.assertEqual(r,1);self.assertEqual(c,1)
  r,c=noise(np.array([10,0,9,100.]));self.assertLess(r,0);self.assertTrue(np.isnan(c))
 def test_effect_bins_use_screen_not_truth(self):
  d=self.data(groups=1,n=2);d.loc[:,['audit_screen_A','audit_screen_B','audit_screen_C']]=np.repeat(np.array([1.2,.7])[:,None],3,axis=1)
  d.loc[:,['audit_eval_A','audit_eval_B','audit_eval_C']]=-1
  a=errors(d,'enformer','strong');b=errors(d,'enformer','weak')
  self.assertTrue(np.isfinite(a[0]));self.assertTrue(np.isnan(a[1]));self.assertTrue(np.isnan(b[0]));self.assertTrue(np.isfinite(b[1]))
 def test_no_label_based_selection(self):
  d=self.data();d['cluster']=d.genomic_block;s=select(d,d.risk_enformer);p=d.copy();p.effect=1000
  for lab in 'ABC':p['audit_screen_'+lab]=1000;p['audit_eval_'+lab]=-1000
  np.testing.assert_array_equal(s,select(p,p.risk_enformer))
 def test_exact_metadata_retention(self):
  a=DatasetAudit(self.data(),'Abell',draws=120);a.run();r=pd.DataFrame(a.retention)
  self.assertTrue(r.equal.all());self.assertTrue((r.policy_retained==np.floor(.25*r.candidate_count)).all())
 def test_unmatched_absent_not_invented(self):
  a=DatasetAudit(self.data(),'Abell',draws=120);a.run();q=pd.DataFrame(a.records)
  h=q.loc[q.dimension.eq('cell_match')&q.level.eq('unmatched')]
  self.assertFalse(h.evaluable.any());self.assertTrue(h.p_raw.eq(1).all());self.assertTrue(h.estimate.isna().all())
 def test_exact70_tests_end_to_end_no_scientific_files(self):
  rows=[]
  for ds in ['Abell','Ulirsch']:
   a=DatasetAudit(self.data(ds),ds,draws=150);rows.extend(a.run())
  self.assertEqual(len(rows),70);p=holm([r['p_raw'] for r in rows]);self.assertTrue(((p>=0)&(p<=1)).all())
 def test_family_isolation_source_lofo(self):
  d=self.data('Kircher',groups=8,n=30);r,m=source_lofo(d,'borzoi');self.assertTrue(np.isfinite(r).all())
  for item in m:self.assertNotIn(item['heldout'],item['fit_groups']+item['calibration_groups'])
 def test_fit_uses_only_source_values(self):
  d=self.data('Kircher',groups=8,n=30);p=fit(d,'borzoi');e=self.data();f=e.copy();f.effect=9999;f.noise_variance=9999
  np.testing.assert_array_equal(predict(e,p),predict(f,p))
 def test_strict_cell_aliases(self):
  import re
  a=json.loads((ROOT/'config/cell_aliases.json').read_text());self.assertFalse(bool(re.search(a['HEK293T'],'DNASE:HEK293',re.I)));self.assertTrue(bool(re.search(a['K562'],'DNASE:K562',re.I)))
 def test_full_cpu_pipeline_temporary_synthetic_IO(self):
  from unittest.mock import patch
  from run_audit import run
  d=pd.concat([self.data('Abell'),self.data('Ulirsch'),self.data('Kircher',groups=8,n=30)],ignore_index=True)
  d['variant_id']=d.row_id;d['condition']='synthetic';d['assay_eligible']=True;d['original_enformer_risk']=d.risk_enformer
  old={}
  for ds in ['Abell','Ulirsch']:
   q=d.loc[d.dataset.eq(ds)];a=DatasetAudit(q,ds,draws=20);r=a.gain_record('enformer','global','all','all')
   old[ds]={'primary_error_reduction':r['estimate'],'retained_audited_variants_candidate':r['policy_audited_calls'],'retained_audited_variants_baseline':r['magnitude_audited_calls']}
  provenance={'mock_only':True,'original_primary_summary':{'decision':'FAIL','datasets':old}}
  with tempfile.TemporaryDirectory() as tmp:
   base=Path(tmp)/'base';base.mkdir();out=Path(tmp)/'out'
   with patch('run_audit.input_reader',return_value=(d,provenance)),patch('run_audit.annotate',side_effect=lambda data,base:(data,pd.DataFrame([{'mock_only':True}]))):run(base,out)
   t=pd.read_csv(out/'subgroup_tests.csv');self.assertEqual(int(t.formal_external.sum()),70)
   self.assertTrue((out/'figure_noise_ceiling.pdf').exists());self.assertTrue((out/'original_primary_reconciliation.csv').exists())
   summary=json.loads((out/'follow_up_summary.json').read_text());self.assertFalse(summary['primary_verdict_changed'])
   # The TemporaryDirectory removes every mock result before the test returns.
 def test_runtime_has_no_Torch_or_API(self):
  self.assertNotIn('torch',sys.modules);self.assertFalse(any(n.startswith('alphagenome') for n in sys.modules))
if __name__=='__main__':unittest.main()
