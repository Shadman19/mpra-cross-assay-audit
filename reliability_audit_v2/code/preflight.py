"""Fail early on missing contexts, incompatible ports or output-crop mismatch.

Metadata checks run without a GPU. Native GPU tests run once per model/code/
weights/hardware identity; they are not efficacy results. No model label is used.
"""
import argparse,json,os,time
os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG',':4096:8')
import numpy as np
import pandas as pd
from common import ROOT,SEED,sha,atomic_json,code_hashes,verify_sources,load_protocol
from models import metadata,select_tracks,MODELS,prepare_weights,load_model,infer_profile,one_hot

def metadata_check():
    verify_sources();load_protocol();d=pd.read_csv(ROOT/'data/variants.csv')
    rows=[]
    for model in MODELS:
        m=metadata(model)
        for e,c in d[['element','condition']].drop_duplicates().itertuples(index=False):
            for mod in ['DNASE','ATAC','CAGE','H3K27ac']:
                ids,rank=select_tracks(m,e,c,mod)
                rows.append({'model':model,'element':e,'condition':c,'modality':mod,'tracks':len(ids),'priority':rank,'track_indices':';'.join(map(str,ids))})
                if mod=='DNASE' and not len(ids):raise ValueError('No predeclared DNASE tracks: '+model+' / '+c)
    pd.DataFrame(rows).to_csv(ROOT/'results/track_metadata_preflight.csv',index=False)
    print('PASS: metadata-only tracks and public raw checksums verified')

def gpu_check(model):
    import torch,py2bit
    if not torch.cuda.is_available():raise RuntimeError('Native GPU preflight needs Wahab gpu partition')
    torch.manual_seed(SEED);torch.set_num_threads(4);torch.use_deterministic_algorithms(True)
    torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False
    torch.backends.cudnn.benchmark=False;torch.backends.cudnn.deterministic=True
    identity={'model_spec':MODELS[model],'code_sha256':code_hashes(),'torch':torch.__version__,'device':torch.cuda.get_device_name(),'protocol':sha(ROOT/'config/protocol.json')}
    target=ROOT/'runtime'/('preflight_'+model+'.json')
    if target.exists() and json.loads(target.read_text()).get('identity')==identity:
        print('Native GPU contract already verified for this identity');return
    folder=prepare_weights(model);net=load_model(model,folder,'cuda')
    ref=py2bit.open(str(ROOT/'runtime/hg38.2bit'));d=pd.read_csv(ROOT/'data/variants.csv');r=d.loc[d.dataset.eq('Kircher') & d.primary].iloc[0]
    center=int(r.pos)-1;length=MODELS[model]['length'];seq=ref.sequence(r.chrom,center-length//2,center+length//2).upper();ref.close()
    ids=select_tracks(metadata(model),r.element,r.condition,'DNASE')[0]
    a=infer_profile(net,model,seq,'cuda',ids);b=infer_profile(net,model,seq,'cuda',ids)
    np.testing.assert_array_equal(a,b) # Reproducibility check, actual weights.
    crop_error=0.0
    if model=='borzoi':
        # Retain the full backbone embedding once; compare both pointwise heads
        # on it, avoiding an expensive second backbone or giant full track output.
        x=torch.from_numpy(one_hot(seq)).T.unsqueeze(0).to('cuda')
        net.crop.target_length=6144
        with torch.inference_mode():
            em=net.get_embs_after_crop(x)
            full=net.final_joined_convs(em)
            full=net.final_softplus(net.human_head(full)[:,ids,:])
            mid=(full.shape[-1]-512)//2
            small=net.final_joined_convs(em[:,:,mid:mid+512])
            small=net.final_softplus(net.human_head(small)[:,ids,:])
            crop_error=float((full[:,:,mid:mid+512]-small).abs().max().cpu())
            torch.testing.assert_close(full[:,:,mid:mid+512],small,rtol=1e-5,atol=1e-5)
        net.crop.target_length=512
        # Same center result through the public full forward method.
        c=infer_profile(net,model,seq,'cuda',ids)
        np.testing.assert_allclose(c,a,rtol=1e-5,atol=1e-5)
    atomic_json({'identity':identity,'native_GPU_verified':True,'deterministic_repeat':True,'Borzoi_pointwise_crop_max_abs_error':crop_error,'output_shape':list(a.shape),'GPU_memory_total_bytes':torch.cuda.get_device_properties(0).total_memory,'peak_allocated_bytes':torch.cuda.max_memory_allocated(),'scientific_efficacy_test':False},target)
    print('PASS: native GPU interface, deterministic repeat and center-crop contract',model)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--metadata-only',action='store_true');p.add_argument('--gpu',action='store_true');p.add_argument('--model',choices=list(MODELS));a=p.parse_args()
    if a.gpu:
        if not a.model:p.error('--gpu needs --model')
        gpu_check(a.model)
    else:metadata_check()
