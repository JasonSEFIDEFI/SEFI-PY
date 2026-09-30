"""Reproduce grid/domain refinements, preserving every original checkpoint."""
import argparse
from pathlib import Path
import subprocess
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output',type=Path,default=Path('broad_runs'))
    ap.add_argument('--data',type=Path,default=Path(__file__).resolve().parents[2]/'sim/vortex_console/data')
    ap.add_argument('--jobs',nargs='+',default=['L16_h040','L24_h050','L24_h040','L32_h050'])
    args=ap.parse_args()
    jobs={'L16_h040':(40,.4,'charged_ring_refined'),
          'L24_h050':(48,.5,'charged_ring_refined'),
          'L24_h040':(60,.4,'charged_ring_refined'),
          'L32_h050':(64,.5,'charged_ring_wide')}
    for name in args.jobs:
        n,h,source=jobs[name]
        subprocess.run([sys.executable,str(Path(__file__).with_name('broad_scan.py')),
            str(args.data/(source+'.npz')),'--output',str(args.output/name),
            '--n',str(n),'--h',str(h),'--m','0','1','2','3','4','--k','12'],check=True)


if __name__=='__main__':
    main()
