"""Apply the scanner's vector retention policy to early campaign outputs."""
import argparse
import json
from pathlib import Path
import numpy as np


def compact(root):
    for path in root.rglob('m*_*.npz'):
        with np.load(path) as data:
            if 'saved_indices' in data:
                continue
            ev,vec=data['eigenvalues'],data['eigenvectors']
        report=json.loads(path.with_suffix('.json').read_text())
        first=path.stem.endswith('_00')
        keep=[j for j,r in enumerate(report['eigenvalues']) if first or
            (r['candidate'] and max(r['weighted_symmetry_overlaps'].values(),default=0)<.9)]
        np.savez_compressed(path,eigenvalues=ev,saved_indices=np.array(keep,dtype=int),eigenvectors=vec[:,keep])


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('root',type=Path)
    compact(ap.parse_args().root)
