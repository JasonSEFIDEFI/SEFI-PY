"""Reproduce the calibration suite from retained profiles (sequential jobs)."""
import argparse
from pathlib import Path
import subprocess
import sys

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',type=Path,required=True)
    data=Path(__file__).resolve().parents[2]/'sim/vortex_console/data'
    if not data.exists():data=Path(__file__).with_name('reference_profiles')
    ap.add_argument('--data',type=Path,default=data)
    a=ap.parse_args()
    solver=Path(__file__).with_name('full_z_spectral.py')
    jobs=[('L16_h050','charged_ring_refined',32,.5,16,.0001),
          ('L16_h040','charged_ring_refined',40,.4,16,.0001),
          ('L24_h0667','charged_ring_wide',36,2/3,16,.0001),
          ('L24_h050','charged_ring_refined',48,.5,16,.0001),
          ('L24_h050_shift001','charged_ring_refined',48,.5,16,.01),
          ('L24_h040','charged_ring_refined',60,.4,12,.01),
          ('L32_h050','charged_ring_wide',64,.5,12,.0001)]
    for name,source,n,h,k,shift in jobs:
        subprocess.run([sys.executable,str(solver),str(a.data/(source+'.npz')),
                        '--output',str(a.output/name),'--n',str(n),'--h',str(h),
                        '--k',str(k),'--shift',str(shift)],check=True)

if __name__=='__main__':main()
