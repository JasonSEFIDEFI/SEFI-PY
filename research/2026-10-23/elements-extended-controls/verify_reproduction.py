from pathlib import Path
import csv
import json
import os
import subprocess
import sys
import numpy as np
ROOT=Path(__file__).resolve().parent
paths=['results/mass_predictions.csv','results/lifetime_predictions.csv','results/historical/predictions.csv']

def load(name):
    return list(csv.DictReader((ROOT/name).open()))

def main():
    before={p:load(p) for p in paths}
    env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1'
    for script in ['analyze.py','historical.py']:
        result=subprocess.run([sys.executable,str(ROOT/script)],capture_output=True,text=True,env=env,check=True)
        (ROOT/(script.replace('.py','')+'-reproduction.log')).write_text(result.stdout+result.stderr)
    checks=[]
    for path in paths:
        a=before[path];b=load(path);assert len(a)==len(b)
        maximum=0.
        for old,new in zip(a,b):
            assert old.keys()==new.keys()
            for key in old:
                try:
                    x=float(old[key]);y=float(new[key]);difference=abs(x-y)
                    assert np.isclose(x,y,rtol=1e-12,atol=1e-10),(path,key,x,y)
                    maximum=max(maximum,difference)
                except ValueError:assert old[key]==new[key]
        checks.append(dict(file=path,rows=len(a),max_numeric_difference=maximum,passed=True))
    (ROOT/'results/reproducibility.json').write_text(json.dumps(checks,indent=2))
    print(json.dumps(checks))

if __name__=='__main__':main()
