"""Shared reproducibility helpers; no predictor or outcome is silently imputed."""
from pathlib import Path
import hashlib,json,os,urllib.request,time
ROOT=Path(__file__).resolve().parents[1]
SEED=20261012
MODS=['DNASE','ATAC','CAGE','H3K27ac']
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def atomic_json(obj,path):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    t=path.with_suffix('.partial');t.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');t.replace(path)
def load_protocol(root=ROOT):
    p=root/'config/protocol.json';lock=json.loads((root/'config/protocol_lock.json').read_text())
    if sha(p)!=lock['sha256']:raise ValueError('Locked protocol changed. Use a new version, never alter this run.')
    return json.loads(p.read_text())
def code_hashes():return {p.name:sha(p) for p in sorted(Path(__file__).parent.glob('*.py'))}
def fetch_checked(url,path,expected,algorithm='sha256'):
    """Atomic download; never accept an HTML challenge page as a data file."""
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():
        h=hashlib.new(algorithm)
        with path.open('rb') as f:
            while b:=f.read(8*1024*1024):h.update(b)
        if h.hexdigest()==expected:return path
        raise ValueError('Cached resource checksum mismatch: '+path.name)
    temp=path.with_suffix('.partial');h=hashlib.new(algorithm)
    for attempt in range(3):
        try:
            h=hashlib.new(algorithm)
            with urllib.request.urlopen(url,timeout=90) as r,temp.open('wb') as f:
                while b:=r.read(8*1024*1024):f.write(b);h.update(b)
            if h.hexdigest()!=expected:raise ValueError('Downloaded checksum mismatch: '+path.name)
            temp.replace(path);return path
        except Exception:
            if attempt==2:raise
            time.sleep(2**attempt)
def verify_sources(root=ROOT):
    for r in json.loads((root/'data/source_manifest.json').read_text()):
        fetch_checked(r['url'],root/r['path'],r['sha256'])
def permission_for_alpha():
    if os.environ.get('ALPHAGENOME_CALIBRATION_PERMISSION')!='1':
        raise PermissionError('Any fitted analysis using AlphaGenome outputs requires applicable output-use permission. Do not set ALPHAGENOME_CALIBRATION_PERMISSION=1 merely to bypass this gate.')
