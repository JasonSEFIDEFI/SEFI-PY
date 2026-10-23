"""Download only versioned public evaluated nuclear data; pin hashes and sources."""
from pathlib import Path
import datetime
import hashlib
import json
import urllib.request

ROOT=Path(__file__).resolve().parent
SOURCES={
 'mass16.txt':'https://www-nds.iaea.org/amdc/ame2016/mass16.txt',
 'mass_1.mas20.txt':'https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt',
 'nubase_4.mas20.txt':'https://www-nds.iaea.org/amdc/ame2020/nubase_4.mas20.txt'}

def main():
    folder=ROOT/'data';folder.mkdir(exist_ok=True)
    meta=[]
    for name,url in SOURCES.items():
        path=folder/name
        if not path.exists():
            request=urllib.request.Request(url,headers={'User-Agent':'SEFI research data audit; Python urllib'})
            with urllib.request.urlopen(request,timeout=90) as response:data=response.read()
            if b'<html' in data[:200].lower():raise ValueError('Unexpected HTML response')
            path.write_bytes(data)
        data=path.read_bytes()
        meta.append(dict(file=name,url=url,sha256=hashlib.sha256(data).hexdigest(),bytes=len(data),
                         accessed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        print(json.dumps(meta[-1]),flush=True)
    (folder/'provenance.json').write_text(json.dumps(meta,indent=2))

if __name__=='__main__':main()
