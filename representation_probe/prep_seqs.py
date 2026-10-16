"""Build REF and ALT insert sequences once, as 0-3 base indices, for target-assay fine-tuning.
Each insert is a variant-centred window carrying the published MPRA flanks out to the model length.
Usage: python prep_seqs.py --genome hg38.2bit --rows CSV.GZ --boda D --dataset Abell --insert-len 200 --model-len 600 --out NPZ"""
import argparse, importlib, sys, types
from pathlib import Path
import numpy as np, pandas as pd, py2bit, torch

ap = argparse.ArgumentParser()
for k in ['genome', 'rows', 'boda', 'out']: ap.add_argument('--' + k, required=True)
ap.add_argument('--dataset', default='Abell'); ap.add_argument('--insert-len', dest='ins', type=int, default=200)
ap.add_argument('--model-len', dest='mlen', type=int, default=600)
a = ap.parse_args()
lt = types.ModuleType('lightning'); pl = types.ModuleType('lightning.pytorch')
pl.LightningModule = torch.nn.Module
def _stub(name):
    if name.startswith('__'): raise AttributeError(name)
    return torch.nn.Module
pl.__getattr__ = _stub; lt.pytorch = pl
sys.modules.setdefault('lightning', lt); sys.modules.setdefault('lightning.pytorch', pl)
BP = Path(a.boda) / 'boda'
for nm, sub in [('bd', ''), ('bd.common', 'common'), ('bd.model', 'model')]:
    m = types.ModuleType(nm); m.__path__ = [str(BP / sub)]; sys.modules.setdefault(nm, m)
C = importlib.import_module('bd.common.constants')
lp = (a.mlen - a.ins) // 2; rp = a.mlen - a.ins - lp
FL, FR = C.MPRA_UPSTREAM[-lp:], C.MPRA_DOWNSTREAM[:rp]
IDX = {c: i for i, c in enumerate('ACGT')}
def enc(s): return np.array([IDX.get(c, 255) for c in s[:a.mlen]], np.uint8)

d = pd.read_csv(a.rows, low_memory=False)
d = d[d.primary.astype(str).str.lower().eq('true') & d.dataset.eq(a.dataset)].copy()
g = py2bit.open(a.genome); R, A, ids, bad = [], [], [], 0
for r in d.itertuples():
    p0 = int(r.pos) - 1; s0 = p0 - 150
    w = g.sequence(str(r.chrom), s0, s0 + 500).upper()
    rf, al = str(r.ref).upper(), str(r.alt).upper(); o = p0 - s0
    if w[o:o + len(rf)] != rf: bad += 1
    wa = w[:o] + al + w[o + len(rf):]
    st = 150 - a.ins // 2; c0 = st + (len(al) - len(rf)) // 2
    R.append(enc(FL + w[st:st + a.ins] + FR)); A.append(enc(FL + wa[c0:c0 + a.ins] + FR)); ids.append(str(r.row_id))
g.close()
np.savez_compressed(a.out, ref=np.array(R), alt=np.array(A), row_id=np.array(ids, dtype=str))
print('wrote', a.out, len(R), 'pairs, ref mismatches', bad, 'len', len(R[0]))
