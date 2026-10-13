"""Resumable GPU scoring. Acquires every mapped source contrast.

Two models run as separate Slurm array tasks, never concurrently in GPU memory.
Fixed2kb genomic tiles reuse reference passes within saturation regions. Every
ALT receives a real model pass; no embedding/linear surrogate is substituted.
"""
import argparse,json,time,hashlib,signal,os
os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
from pathlib import Path
import numpy as np
import pandas as pd
import py2bit
from common import ROOT,SEED,MODS,load_protocol,sha,atomic_json,code_hashes
from models import metadata,select_tracks,MODELS,prepare_weights,load_model,infer_profile,reverse_complement,strand_permutation,centered_scores
STOP=False
def stop_handler(*args):
    global STOP;STOP=True
def genome_window(ref,chrom,start,width):
    size=ref.chroms(chrom)
    a=max(0,start);b=min(size,start+width)
    mid=ref.sequence(chrom,a,b).upper() if b>a else ''
    left=min(width,max(0,-start))
    seq='N'*left+mid+'N'*(width-left-len(mid))
    if len(seq)!=width:raise ValueError('Reference window length mismatch')
    return seq
def edit_contrast(genomic,offset,ref,alt,tail):
    # Explicit library REF contrasts are preserved, even when source REF is not
    # hg38. A deletion consumes one base and adds right genomic context.
    reference=(genomic[:offset]+ref+genomic[offset+len(ref):])[:len(genomic)]
    alternate=reference[:offset]+('' if alt=='-' else alt)+reference[offset+len(ref):]
    alternate=(alternate+tail)[:len(genomic)]
    if len(reference)!=len(genomic) or len(alternate)!=len(genomic):raise ValueError('Edited input length mismatch')
    return reference,alternate
def run(model,shard=0,shards=8,max_seconds=19800,root=ROOT):
    import torch
    protocol=load_protocol(root)
    if not torch.cuda.is_available():raise RuntimeError('This job must run on Wahab gpu partition; CUDA unavailable')
    torch.manual_seed(SEED);np.random.seed(SEED);torch.set_num_threads(4)
    torch.backends.cudnn.benchmark=False;torch.backends.cudnn.deterministic=True
    torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False
    torch.use_deterministic_algorithms(True,warn_only=False)
    signal.signal(signal.SIGTERM,stop_handler);signal.signal(signal.SIGUSR1,stop_handler)
    d=pd.read_csv(root/'data/variants.csv')
    # Stable hashof context groups keeps neighboring variants/reference passes
    # together and gives the same partition on every resume.
    d=d.loc[d.context_id.map(lambda s:int(hashlib.sha256(s.encode()).hexdigest()[:12],16)%shards).eq(shard)].copy()
    out=root/'scores'/model;out.mkdir(parents=True,exist_ok=True)
    identity={'model':model,'model_spec':MODELS[model],'protocol_sha256':sha(root/'config/protocol.json'),'source_rows_sha256':sha(root/'data/variants.csv'),'code_sha256':code_hashes(),'shard':shard,'shards':shards,'dtype':'float32','torch_version':torch.__version__,'cuda_version':torch.version.cuda,'cudnn_version':torch.backends.cudnn.version()}
    state=out/f'shard{shard:02d}.json';journal=out/f'shard{shard:02d}.jsonl'
    if state.exists() and json.loads(state.read_text())['identity']!=identity:raise ValueError('Checkpoint identity changed; keep a separate run version')
    completed={}
    if journal.exists():
        for line in journal.read_text().splitlines():
            # A terminated write can leave only the final line partial.
            try:r=json.loads(line)
            except json.JSONDecodeError:break
            if r['row_id'] in completed:raise ValueError('Duplicate checkpoint row')
            completed[r['row_id']]=r
        # Rewrite valid entries atomically before appending after a partial line.
        tmp=journal.with_suffix('.partial');tmp.write_text(''.join(json.dumps(r,allow_nan=False)+'\n' for r in completed.values()));tmp.replace(journal)
    atomic_json({'identity':identity,'complete':False,'scored_rows':len(completed)},state)
    if set(completed)==set(d.row_id):
        atomic_json({'identity':identity,'complete':True,'scored_rows':len(completed),'real_model_scores':True},state);return
    started=time.monotonic();net=load_model(model,prepare_weights(model,root),'cuda')
    meta=metadata(model,root);pair=strand_permutation(meta);spec=MODELS[model]
    # Allselected tracks and their strand partners are fixed from metadata.
    selections={};selectionrows=[];union=set()
    for (element,condition),g in d.groupby(['element','condition']):
        selected={}
        for mod in MODS:
            ids,rank=select_tracks(meta,element,condition,mod);selected[mod]=ids
            union.update(ids);union.update(pair[ids])
            selectionrows.append({'element':element,'condition':condition,'modality':mod,'n_tracks':len(ids),'proxy_priority':rank,'track_indices':';'.join(map(str,ids)),'track_identifiers':';'.join(meta.iloc[ids].identifier.astype(str))})
        selections[(element,condition)]=selected
        # Source-quality-independent missingness cannot silently drop a locus.
        if not len(selected['DNASE']):raise ValueError('No declared DNASE context for '+condition)
    indices=np.array(sorted(union),int);lookup={v:i for i,v in enumerate(indices)}
    rcorder=np.array([lookup[pair[i]] for i in indices])
    pd.DataFrame(selectionrows).to_csv(out/f'tracks{shard:02d}.csv',index=False)
    ref=py2bit.open(str(root/'runtime/hg38.2bit'))
    # A journal row is committed immediately after each condition's summaries.
    # Full profile arrays never enter the ZIP or become substitute predictions.
    with journal.open('a') as log:
        for context,g in d.sort_values(['context_id','variant_id','row_id']).groupby('context_id',sort=True):
            todo=g.loc[~g.row_id.isin(completed)]
            if not len(todo):continue
            if STOP or time.monotonic()-started>max_seconds:break
            tile=int(g.tile_start.iloc[0])+1024;start=tile-spec['length']//2
            genomic=genome_window(ref,g.chrom.iloc[0],start,spec['length'])
            # Referenceonly two passes per tile. Library REF mismatches get
            # their own passes below instead of changing the provided allele.
            rf=infer_profile(net,model,genomic,'cuda',indices)
            rc=infer_profile(net,model,reverse_complement(genomic),'cuda',indices)[::-1][:,rcorder]
            for vid,vrows in todo.groupby('variant_id',sort=True):
                if STOP or time.monotonic()-started>max_seconds:break
                r=vrows.iloc[0];offset=int(r.pos)-1-start
                tail=genome_window(ref,r.chrom,start+spec['length'],max(1,len(r.ref)))
                rs,als=edit_contrast(genomic,offset,r.ref,r.alt,tail)
                if rs==genomic:rforward,rreverse=rf,rc
                else:
                    rforward=infer_profile(net,model,rs,'cuda',indices)
                    rreverse=infer_profile(net,model,reverse_complement(rs),'cuda',indices)[::-1][:,rcorder]
                af=infer_profile(net,model,als,'cuda',indices)
                ar=infer_profile(net,model,reverse_complement(als),'cuda',indices)[::-1][:,rcorder]
                scores=centered_scores(rforward,af,rreverse,ar,int(r.pos)-1-tile,spec['bin_bp'],meta,indices)
                for row in vrows.itertuples():
                    result={'row_id':row.row_id,'variant_id':row.variant_id,'model':model,'real_model_scores':True,'scoring_device':torch.cuda.get_device_name(),'cuda_version':torch.version.cuda}
                    for mod,ids in selections[(row.element,row.condition)].items():
                        s=scores[[lookup[i] for i in ids]]
                        result[mod+'_mean']=float(s.mean()) if len(s) else None
                        result[mod+'_std']=float(s.std()) if len(s) else None
                        result[mod+'_n']=len(s)
                    log.write(json.dumps(result,allow_nan=False)+'\n');log.flush();completed[row.row_id]=result
            print(model,'shard',shard,'rows',len(completed),'of',len(d),flush=True)
    ref.close()
    complete=set(completed)==set(d.row_id)
    atomic_json({'identity':identity,'complete':complete,'scored_rows':len(completed),'expected_rows':len(d),'real_model_scores':True,'device':torch.cuda.get_device_name(),'wall_seconds_this_resume':time.monotonic()-started,'stop_note':'Resume the exact same shard; do not delete valid checkpoints'},state)
    if not complete:raise SystemExit(75)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--model',choices=['enformer','borzoi'],required=True);p.add_argument('--shard',type=int,default=0);p.add_argument('--shards',type=int,default=8);p.add_argument('--max-seconds',type=int,default=19800);a=p.parse_args();run(a.model,a.shard,a.shards,a.max_seconds)
