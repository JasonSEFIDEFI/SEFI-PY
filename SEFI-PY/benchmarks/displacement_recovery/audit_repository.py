"""Inspect every repository-owned source in the pinned inventory, read-only."""
import ast
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[2]
inventory=json.loads((ROOT/'repository_inventory.json').read_text())
rows=[]; syntax_issues=[]; missing=[]
for entry in inventory['entries']:
    path=entry['path']
    if entry['type']!='blob' or Path(path).suffix not in {'.py','.js','.ts','.tsx','.jsx','.qs','.toml','.ini'}:continue
    local=REPO/path
    if not local.is_file():missing.append(path);continue
    data=local.read_bytes()
    blob=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
    lf=data.replace(b'\r\n',b'\n')
    lf_blob=hashlib.sha1(b'blob '+str(len(lf)).encode()+b'\0'+lf).hexdigest()
    text=data.decode('utf-8-sig',errors='replace')
    status='not Python'; names=[]
    if path.endswith('.py'):
        try:
            tree=ast.parse(text);status='parses'
            names=[n.name for n in ast.walk(tree) if isinstance(n,(ast.FunctionDef,ast.ClassDef))]
        except SyntaxError as err:
            status=f'SyntaxError line {err.lineno}: {err.msg}';syntax_issues.append({'path':path,'error':status})
    rows.append({'path':path,'lines':len(text.splitlines()),'bytes':len(data),
                 'matches_pinned_blob':blob==entry['sha'] or lf_blob==entry['sha'],
                 'match_after_crlf_normalization':blob!=entry['sha'] and lf_blob==entry['sha'],
                 'syntax':status,'definitions':'; '.join(names)})
ROOT.joinpath('audit').mkdir(exist_ok=True)
with (ROOT/'audit/source_inventory.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
result={'pinned_commit':inventory['sha'],'remote_tree_entries':inventory['entry_count'],
        'read_source_files':len(rows),'line_count':sum(x['lines'] for x in rows),
        'matching_source_files':sum(x['matches_pinned_blob'] for x in rows),
        'syntax_errors':syntax_issues,'missing_local_sources':missing,
        'scope':'Every owned source file in recursively enumerated tree; dependencies, caches, videos excluded. Static syntax parsing is not runtime verification. Binary assets enumerated, not scientifically validated.'}
(ROOT/'audit/audit_record.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
