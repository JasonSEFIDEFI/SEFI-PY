"""Run reproducible dense/sparse comparisons; no speed assertions."""
import os
os.environ['OPENBLAS_NUM_THREADS']='2'
os.environ['OMP_NUM_THREADS']='2'
import json,platform,statistics,time
from pathlib import Path
import numpy as np
import scipy
import physics

def main():
    reference,_=physics.load_profile(Path(__file__).parent/'data/charged_ring_refined.npz')
    results=[]
    for n in (32,):
        p=physics.validate(dict(n=n,c=.58))
        samples={k:[] for k in ('dense','sparse')};fields={};reports={}
        for repeat in range(3):
            for backend in ('dense','sparse') if repeat%2==0 else ('sparse','dense'):
                started=time.perf_counter()
                fields[backend],reports[backend]=physics.solve(p,reference,backend=backend)
                samples[backend].append(time.perf_counter()-started)
        results.append(dict(n=n,parameters=p,samples_seconds=samples,
            median_seconds={k:statistics.median(v) for k,v in samples.items()},
            speedup=statistics.median(samples['dense'])/statistics.median(samples['sparse']),
            max_field_difference=max(float(np.max(abs(fields['dense'][k]-fields['sparse'][k]))) for k in ('u','w','s')),
            reports=reports))
    print(json.dumps(dict(python=platform.python_version(),platform=platform.platform(),
        numpy=np.__version__,scipy=scipy.__version__,threads=2,results=results),indent=2))
if __name__=='__main__':main()
