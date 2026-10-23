from pathlib import Path
import hashlib
import json
import zipfile
ROOT=Path(__file__).resolve().parent

def main():
    manifest={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.rglob('*')
        if p.is_file() and '__pycache__' not in p.parts and '.git' not in p.parts and p.name!='SHA256SUMS.json'}
    (ROOT/'SHA256SUMS.json').write_text(json.dumps(manifest,indent=2))
    archive=ROOT.parent/'SEFI_ELEMENTS_Research_Package_2026-10-23.zip'
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for name in sorted([*manifest,'SHA256SUMS.json']):
            z.write(ROOT/name,arcname=ROOT.name+'/'+name)
    with zipfile.ZipFile(archive) as z:assert z.testzip() is None
    print(json.dumps(dict(archive=str(archive),files=len(manifest)+1,bytes=archive.stat().st_size)))

if __name__=='__main__':main()
