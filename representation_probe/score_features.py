"""GPU feature extraction for the preregistered representation probe (Enformer only).

For every primary audit row: B = full 5,313-track ALT-REF log2 score; C = 3,072-d trunk embedding
ALT-REF averaged over the four centre bins; D = 6-mer count difference in a 21 bp window.
Forward and reverse-complement passes are averaged. Reuses the audit's pinned model, genome and
sequence-editing code unchanged. Resumable per shard (one .npz per finished context group).
Usage: python score_features.py --base <Reliability_Audit_v2 dir> --rows <annotated_audit_variants.csv.gz>
       --out <dir> --shard I --shards N
"""
import argparse, hashlib, os, sys, time, signal
from itertools import product
from pathlib import Path
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('--base', required=True); ap.add_argument('--rows', required=True); ap.add_argument('--out', required=True)
ap.add_argument('--shard', type=int, default=0); ap.add_argument('--shards', type=int, default=8)
ap.add_argument('--max-seconds', type=int, default=19200)
a = ap.parse_args()
BASE = Path(a.base).resolve()
sys.path.insert(0, str(BASE / 'code'))
os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG', ':4096:8')
import torch, py2bit  # noqa: E402
from models import MODELS, prepare_weights, load_model, one_hot, reverse_complement, metadata, centered_scores  # noqa: E402
from score_open_models import genome_window, edit_contrast  # noqa: E402

STOP = False
def _stop(*_):
    global STOP; STOP = True
signal.signal(signal.SIGTERM, _stop); signal.signal(signal.SIGUSR1, _stop)

KMERS = {''.join(p): i for i, p in enumerate(product('ACGT', repeat=6))}
def kmer_vec(s):
    v = np.zeros(4096, np.float32)
    for i in range(len(s) - 5):
        j = KMERS.get(s[i:i + 6])
        if j is not None: v[j] += 1
    return v

def run_pass(net, seq):
    x = torch.from_numpy(one_hot(seq)).unsqueeze(0).to('cuda')
    with torch.inference_mode():
        out, emb = net(x, head='human', return_embeddings=True)
    out = out[0].float().cpu().numpy(); emb = emb[0].float().cpu().numpy()
    if out.shape[-1] != 5313 or emb.shape[-1] != 3072: raise ValueError(f'Unexpected shapes {out.shape} {emb.shape}')
    if not (np.isfinite(out).all() and np.isfinite(emb).all()): raise ValueError('Nonfinite output')
    return out, emb

def both(net, seq):
    f_out, f_emb = run_pass(net, seq)
    r_out, r_emb = run_pass(net, reverse_complement(seq))
    return f_out, f_emb, r_out[::-1], r_emb[::-1]   # Enformer human tracks are unstranded

def main():
    if not torch.cuda.is_available(): raise RuntimeError('Run on the Wahab gpu partition')
    torch.manual_seed(20261012); torch.backends.cudnn.deterministic = True; torch.backends.cudnn.benchmark = False
    torch.backends.cuda.matmul.allow_tf32 = False; torch.backends.cudnn.allow_tf32 = False
    d = pd.read_csv(a.rows, low_memory=False)
    d = d[d.primary.astype(str).str.lower().eq('true')].copy()
    d = d[d.context_id.map(lambda s: int(hashlib.sha256(str(s).encode()).hexdigest()[:12], 16) % a.shards).eq(a.shard)]
    out = Path(a.out) / f'shard{a.shard:02d}'; out.mkdir(parents=True, exist_ok=True)
    spec = MODELS['enformer']; L = spec['length']; bp = spec['bin_bp']
    net = load_model('enformer', prepare_weights('enformer', BASE), 'cuda')
    meta = metadata('enformer', BASE); allidx = np.arange(spec['tracks'])
    ref = py2bit.open(str(BASE / 'runtime/hg38.2bit')); t0 = time.monotonic(); done = 0
    for ctx, g in d.sort_values(['context_id', 'row_id']).groupby('context_id', sort=True):
        f = out / (hashlib.sha256(str(ctx).encode()).hexdigest()[:16] + '.npz')
        if f.exists(): done += len(g); continue
        if STOP or time.monotonic() - t0 > a.max_seconds: break
        tile = int(g.tile_start.iloc[0]) + 1024; start = tile - L // 2
        genomic = genome_window(ref, g.chrom.iloc[0], start, L)
        base_pass = both(net, genomic)
        rows, B, C, D = [], [], [], []
        for vid, vr in g.groupby('variant_id', sort=True):
            r = vr.iloc[0]; offset = int(r.pos) - 1 - start
            tail = genome_window(ref, r.chrom, start + L, max(1, len(str(r.ref))))
            rs, als = edit_contrast(genomic, offset, str(r.ref), str(r.alt), tail)
            rp = base_pass if rs == genomic else both(net, rs)
            ap_ = both(net, als)
            b = centered_scores(rp[0], ap_[0], rp[2], ap_[2], int(r.pos) - 1 - tile, bp, meta, allidx)
            n = rp[1].shape[0]; mid = n // 2 + int((int(r.pos) - 1 - tile) // bp); h = 256 // bp
            c = ((ap_[1][mid - h:mid + h] - rp[1][mid - h:mid + h]).mean(0) + (ap_[3][mid - h:mid + h] - rp[3][mid - h:mid + h]).mean(0)) / 2
            dk = kmer_vec(als[offset - 10:offset + 11]) - kmer_vec(rs[offset - 10:offset + 11])
            for row in vr.itertuples():
                rows.append(row.row_id); B.append(b.astype(np.float16)); C.append(c.astype(np.float16)); D.append(dk.astype(np.int8))
        np.savez(f.with_suffix('.tmp.npz'), row_id=np.array(rows, dtype=str), B=np.array(B), C=np.array(C), D=np.array(D))
        os.replace(f.with_suffix('.tmp.npz'), f); done += len(g)
        print('shard', a.shard, 'rows', done, 'of', len(d), flush=True)
    ref.close()
    if done < len(d): raise SystemExit(75)
    (out / 'COMPLETE').write_text(str(done))

if __name__ == '__main__':
    main()
