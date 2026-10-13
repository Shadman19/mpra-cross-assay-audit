"""Pinned released open-model ports and explicit track/strand contracts."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from common import ROOT,MODS,fetch_checked
import sys
sys.path.insert(0,str(ROOT/'vendor'))
MODELS={
 'enformer':{'repo':'EleutherAI/enformer-official-rough','revision':'affe5713ae9017460706a44108289b13c5fee16c','weights':'pytorch_model.bin','sha256':'99b09d602e195d89c7d4debe144bb2f43907ba0d74006e97098e99d9171c439c','length':196608,'bin_bp':128,'tracks':5313},
 'borzoi':{'repo':'johahi/borzoi-replicate-0','revision':'521db103372769581b8d1d616e6a3b496da490c5','weights':'model.safetensors','sha256':'d839593c3350bc3aea93b7a6d1fc3ce2b1629e0f4322ecb3c6b6743d14932f75','length':524288,'bin_bp':32,'tracks':7611}}
PATTERNS={
 'BCL11A':[r'\bHEL\b',r'K-?562',r'erythro'], 'F9':[r'Hep-?G2',r'liver|hepatocyte'],
 'FOXE1':[r'HeLa',r'cervi'], 'GP1BB':[r'\bHEL\b',r'K-?562',r'erythro'],
 'HBB':[r'\bHEL\b',r'K-?562',r'erythro'],'HBG1':[r'\bHEL\b',r'K-?562',r'erythro'],
 'HNF4A':[r'HEK-?293T',r'HEK-?293'],'IRF4':[r'SK[- ]?MEL[- ]?28',r'SK[- ]?MEL',r'melanocyte'],
 'IRF6':[r'HaCaT',r'keratinocyte'],'LDLR':[r'Hep-?G2',r'liver|hepatocyte'],
 'MSMB':[r'HEK-?293T',r'HEK-?293'],'MYCrs11986220':[r'LNCaP',r'prostate'],
 'MYCrs6983267':[r'HEK-?293T',r'HEK-?293'],'PKLR':[r'K-?562',r'erythro'],
 'RET':[r'SH[- ]?SY5Y|SK[- ]?N[- ]?SH',r'neural|neuron|brain'],'SORT1':[r'Hep-?G2',r'liver|hepatocyte'],
 'TCF7L2':[r'beta[ _-]?cell|islet',r'pancreas'],'TERT':[r'HEK-?293T',r'HEK-?293'],
 'UC88':[r'SH[- ]?SY5Y|SK[- ]?N[- ]?SH',r'neural|neuron|brain'],'ZFAND3':[r'beta[ _-]?cell|islet',r'pancreas'],
 'ZRS':[r'fibroblast',r'mesenchymal'],'Ulirsch_variant':[r'K-?562'],'Abell_variant':[r'GM-?12878']}
def metadata(model,root=ROOT):
    m=pd.read_csv(root/'data/model_metadata'/(model+'_targets.tsv'),sep='\t')
    if 'Unnamed: 0' in m:m=m.rename(columns={'Unnamed: 0':'index'})
    assert len(m)==MODELS[model]['tracks'] and m['index'].tolist()==list(range(len(m)))
    return m
def select_tracks(meta,element,condition,mod):
    text=meta.description.fillna('')
    mm=text.str.startswith('CHIP:') & text.str.contains('H3K27ac',case=False) if mod=='H3K27ac' else text.str.startswith(mod+':')
    patterns=PATTERNS[element]
    if element=='TERT' and condition in ['TERT-GBM','TERT-GAa','TERT-GSc']:patterns=[r'SF7996',r'glioblastoma|U-?87|U-?251',r'brain|astrocyte']
    for rank,pattern in enumerate(patterns):
        mask=mm & text.str.contains(pattern,case=False,regex=True)
        if mask.any():
            # Do not mix untreated external assays with induced/drug conditions
            # when untreated tracks exist. The metadata rule never sees labels.
            clean=mask & ~text.str.contains('treated|phase|induced|stimulat|knockdown|shRNA',case=False,regex=True)
            if clean.any():mask=clean
            return meta.loc[mask,'index'].to_numpy(int),rank
    return np.array([],int),None
def reverse_complement(seq):return seq.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def one_hot(seq):
    lut=np.zeros((256,4),np.float32)
    for i,n in enumerate(b'ACGT'):lut[n,i]=1
    return lut[np.frombuffer(seq.encode(),dtype=np.uint8)]
def strand_permutation(meta):
    if 'strand_pair' in meta:
        p=meta.strand_pair.to_numpy(int)
        if not np.array_equal(p[p],np.arange(len(p))):raise ValueError('Invalid strand-pair permutation')
        return p
    # Released Enformer human tracks are unstranded. Borzoi CAGE/RNA are not.
    return np.arange(len(meta))
def prepare_weights(model,root=ROOT):
    from huggingface_hub import snapshot_download
    spec=MODELS[model]
    folder=Path(snapshot_download(spec['repo'],revision=spec['revision'],allow_patterns=['config.json',spec['weights']],cache_dir=str(root/'runtime/hf_cache')))
    fetch_checked('https://huggingface.co/'+spec['repo']+'/resolve/'+spec['revision']+'/'+spec['weights'],folder/spec['weights'],spec['sha256'])
    return folder
def load_model(model,folder,device):
    import torch
    if model=='enformer':
        from enformer_pytorch import Enformer
        from enformer_pytorch.modeling_enformer import Attention
        net=Enformer.from_pretrained(str(folder),local_files_only=True)
        # The helper only enables corrected TF gamma for the literal Hub name.
        # Explicitly enable it for a pinned local snapshot as well.
        for layer in net.modules():
            if isinstance(layer,Attention):layer.use_tf_gamma=True
        net.use_tf_gamma=True
    else:
        from borzoi_pytorch import Borzoi
        net=Borzoi.from_pretrained(str(folder),local_files_only=True)
        # Crop only the output embedding before the 1x1 final head. This keeps
        # full524kb context and unchanged center predictions while saving RAM.
        net.crop.target_length=512
    return net.to(device).eval()
def infer_profile(net,model,seq,device,selected=None):
    import torch
    x=torch.from_numpy(one_hot(seq)).unsqueeze(0).to(device)
    if model=='borzoi':x=x.transpose(1,2)
    # Full float32 is the locked scoring precision. No seed-dependent dropout,
    # mixed-precision switching, or changing dtype to rescue a result.
    with torch.inference_mode():
        y=net(x,head='human') if model=='enformer' else net(x,is_human=True)
        if model=='borzoi':y=y.transpose(1,2)
        assert y.ndim==3 and y.shape[-1]==MODELS[model]['tracks']
        if selected is not None:y=y[:,:,selected]
        result=y[0].float().cpu().numpy()
    if not np.isfinite(result).all() or result.min()<0:raise ValueError('Nonfinite/negative track prediction')
    return result
def centered_scores(ref_fwd,alt_fwd,ref_rc,alt_rc,offset,bp,meta,indices):
    """Raw log2 ALT/REF score, 512bp mask, aligned reverse complement bins."""
    n=ref_fwd.shape[0];mid=n//2+int(offset//bp);half=256//bp
    lo=mid-half;hi=mid+half
    if lo<0 or hi>n:raise ValueError('Variant outside center output crop')
    # ref_rc/alt_rc already have bins reversed and stranded targets remapped.
    # Borzoi target squashing is undone pertrack using released metadata.
    def undo(x):
        scale=meta.iloc[indices].scale.to_numpy(float)
        z=x/scale
        if 'clip_soft' in meta:
            # Released non-RNA targets use soft clipping followed by scale.
            # Undo in the reverse order. RNA's separate3/4 power transform is
            # irrelevant here because RNA tracks are not part of this score.
            clip=meta.iloc[indices].clip_soft.to_numpy(float)
            z=np.where(z>clip,(z-clip)**2+clip,z)
        if not np.isfinite(z).all():raise ValueError('Invalid released target scale')
        return z
    r=(undo(ref_fwd[lo:hi]).sum(axis=0)+undo(ref_rc[lo:hi]).sum(axis=0))/2
    a=(undo(alt_fwd[lo:hi]).sum(axis=0)+undo(alt_rc[lo:hi]).sum(axis=0))/2
    return np.log2(a+1)-np.log2(r+1)
