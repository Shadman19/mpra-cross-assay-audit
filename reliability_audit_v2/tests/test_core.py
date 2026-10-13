"""Software tests, including a synthetic smoke audit. No efficacy inference.

All mock predictions and mock results live in a TemporaryDirectory, never in
results/. Native released-weight GPU tests are separately run on Wahab.
"""
import unittest,tempfile,sys,json,shutil,os
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'code'))
from common import load_protocol,verify_sources,permission_for_alpha
from models import metadata,strand_permutation,reverse_complement,centered_scores
from score_open_models import edit_contrast
from prepare_data import paired_effect
from audit import select,deployment_features,fit_policy,apply_policy,paired_bootstrap,noise_corrected,evaluate

class PublicSourceTests(unittest.TestCase):
    def test_checksums_and_protocol(self):verify_sources();load_protocol()
    def test_actual_data_replication_and_independence(self):
        d=pd.read_csv(ROOT/'data/variants.csv');p=d.loc[d.primary];r=pd.read_csv(ROOT/'data/replicates.csv')
        self.assertTrue(d.row_id.is_unique)
        self.assertTrue(p.loc[p.dataset.eq('Abell')].variant_id.is_unique)
        sets={ds:set(p.loc[p.dataset.eq(ds),'genomic_block']) for ds in ['Kircher','Ulirsch','Abell']}
        for a,b in [('Kircher','Ulirsch'),('Kircher','Abell'),('Ulirsch','Abell')]:self.assertFalse(sets[a]&sets[b])
        ab=r.loc[r.component.eq('paired_logratio_effect')].groupby('row_id').effect.agg(['mean','var','size'])
        ae=d.set_index('row_id').loc[ab.index]
        np.testing.assert_allclose(ab['mean'],ae.effect,atol=1e-12);np.testing.assert_allclose(ab['var']/3,ae.noise_variance,atol=1e-12)
        self.assertTrue(ab['size'].eq(3).all())
        ul=r.loc[r.row_id.str.startswith('Ulirsch|')];rna=ul.loc[ul.component.eq('RNA_logratio')].groupby('row_id').effect.agg(['mean','var','size']);dna=ul.loc[ul.component.eq('DNA_logratio')].groupby('row_id').effect.agg(['mean','var','size'])
        ue=d.set_index('row_id').loc[rna.index]
        np.testing.assert_allclose(rna['mean']-dna['mean'],ue.effect,atol=1e-12);np.testing.assert_allclose(rna['var']/6+dna['var']/2,ue.noise_variance,atol=1e-12)
        self.assertTrue(rna['size'].eq(6).all());self.assertTrue(dna['size'].eq(2).all())
    def test_native_imports_and_strand_contract(self):
        sys.path.insert(0,str(ROOT/'vendor'))
        from enformer_pytorch import Enformer
        from borzoi_pytorch import Borzoi
        self.assertTrue(callable(Enformer.from_pretrained));self.assertTrue(callable(Borzoi.from_pretrained))
        for name in ['enformer','borzoi']:
            m=metadata(name);q=strand_permutation(m);np.testing.assert_array_equal(q[q],np.arange(len(q)))

class StatisticalTests(unittest.TestCase):
    @staticmethod
    def fixture():
        rng=np.random.default_rng(2017);rows=[]
        # Synthetic effects are explicitly arbitrary software-test labels.
        for dataset in ['Kircher','Ulirsch','Abell']:
            for group in range(20):
                for j in range(40):
                    e=rng.normal();b=e+rng.normal(0,.15);y=e+rng.normal(0,.4)
                    r={'row_id':dataset+'|'+str(group)+'|'+str(j),'dataset':dataset,'primary':True,'fit_group':'family'+str(group),'genomic_block':dataset+':'+str(group),'effect':y,'noise_variance':.03,'enformer_DNASE_mean':e,'borzoi_DNASE_mean':b,'enformer_DNASE_std':abs(rng.normal(0,.2))}
                    for model in ['enformer','borzoi']:
                        for mod in ['ATAC','CAGE','H3K27ac']:r[model+'_'+mod+'_mean']=e+rng.normal(0,.1)
                    for lab in 'ABC':r['audit_screen_'+lab]=y+rng.normal(0,.1);r['audit_eval_'+lab]=y+rng.normal(0,.1)
                    rows.append(r)
        return pd.DataFrame(rows)
    def test_matched_retention_and_external_label_poison(self):
        d=self.fixture()
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);(root/'config').mkdir();shutil.copy(ROOT/'config/protocol.json',root/'config');shutil.copy(ROOT/'config/protocol_lock.json',root/'config')
            policy=fit_policy(d,root);a=apply_policy(d,policy)
            poison=d.copy();ix=poison.dataset.ne('Kircher');poison.loc[ix,'effect']=10000;poison.loc[ix,'noise_variance']=10000
            for lab in 'ABC':poison.loc[ix,'audit_screen_'+lab]=0;poison.loc[ix,'audit_eval_'+lab]=0
            other=fit_policy(poison,root);b=apply_policy(poison,other)
            np.testing.assert_array_equal(a.risk_raw,b.risk_raw)
            np.testing.assert_array_equal(a.risk_probability,b.risk_probability)
            methods=['Source-trained risk','Magnitude only','Cross-model disagreement','Cross-track disagreement','Random ranking']
            g=a.loc[a.dataset.eq('Abell') & a.genomic_block.eq('Abell:0')]
            for rate in [.1,.25,.5,.75]:self.assertEqual({len(select(g,k,rate)) for k in methods},{int(np.floor(rate*len(g)))})
    def test_paired_block_bootstrap_known_identity(self):
        a=pd.DataFrame({'block':['b1','b2'],'error_sum':[1,2],'audit_calls':[4,4]});v=paired_bootstrap(a,a,np.random.default_rng(3),10000)
        self.assertEqual(v,(0.0,[0.0,0.0],[0.0,0.0]))
    def test_noise_ceiling_no_clipping_or_fake_precision(self):
        g=pd.DataFrame({'effect':np.arange(12),'forecast':np.arange(12),'noise_variance':1000})
        a,rel,corr,status=noise_corrected(g);self.assertLess(rel,0);self.assertTrue(np.isnan(corr))
        g.noise_variance=1;a,rel,corr,status=noise_corrected(g);self.assertGreater(corr,1);self.assertIn('unstable',status)
    def test_alpha_permission_closed(self):
        old=os.environ.pop('ALPHAGENOME_CALIBRATION_PERMISSION',None)
        try:
            with self.assertRaises(PermissionError):permission_for_alpha()
            with self.assertRaises(PermissionError):fit_policy(self.fixture(),alpha=True)
        finally:
            if old is not None:os.environ['ALPHAGENOME_CALIBRATION_PERMISSION']=old
    def test_end_to_end_mock_smoke_outputs_only_temporary(self):
        d=self.fixture()
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);(root/'config').mkdir();shutil.copy(ROOT/'config/protocol.json',root/'config');shutil.copy(ROOT/'config/protocol_lock.json',root/'config')
            policy=fit_policy(d,root);status=evaluate(d,policy,root,boot=80)
            self.assertIn(status['decision'],['PASS','FAIL','INCONCLUSIVE'])
            m=pd.read_csv(root/'results/matched_retention_by_block.csv')
            self.assertTrue(m.groupby(['dataset','block','retention']).retained.nunique().eq(1).all())
            self.assertTrue((root/'results/figure1_external_gain.png').exists())
            self.assertTrue((root/'results/precision_with_bootstrap_CIs.csv').exists())

class SequenceTests(unittest.TestCase):
    def test_real_allele_edit_not_reference_mutation(self):
        ref,alt=edit_contrast('ACGTACGT',2,'G','T','A');self.assertEqual(ref,'ACGTACGT');self.assertEqual(alt,'ACTTACGT')
        ref,alt=edit_contrast('ACGTACGT',2,'G','-','A');self.assertEqual(alt,'ACTACGTA')
        ref,alt=edit_contrast('ACGTACGT',2,'A','T','A');self.assertEqual(ref,'ACATACGT');self.assertEqual(alt,'ACTTACGT')
        self.assertEqual(reverse_complement(reverse_complement('ACGTN')),'ACGTN')
    def test_center_score_and_inverse_track_scale(self):
        m=pd.DataFrame({'scale':[2.,4.]});r=np.ones((8,2))*2;a=r*2
        v=centered_scores(r,a,r,a,0,128,m,np.array([0,1]))
        np.testing.assert_allclose(v,np.log2(np.array([9.,5.]))-np.log2(np.array([5.,3.])))

if __name__=='__main__':unittest.main()
