"""Time-budgeted local research batch. No cloud allocation or console mutation."""
from cylindrical3d import *
from datetime import datetime,timezone

DEADLINE=datetime(2026,10,19,21,36,0,tzinfo=timezone.utc).timestamp()
tasks=[
 ('nonlinear_m2_n16',dict(T=12,ntheta=16,dt=.025,epsilon=.01)),
 ('nonlinear_m2_n24',dict(T=12,ntheta=24,dt=.025,epsilon=.01)),
 ('nonlinear_m2_n16_half_dt',dict(T=12,ntheta=16,dt=.0125,epsilon=.01)),
 ('nonlinear_strong_m2_n16',dict(T=24,ntheta=16,dt=.025,epsilon=.05)),
 ('nonlinear_m2_wide_n16',dict(T=12,ntheta=16,dt=.025,epsilon=.01,reference='wide')),
]
record=[]
for name,kwargs in tasks:
    remaining=DEADLINE-time.time()
    if remaining<45:break
    result=run_case(name,wall_limit=min(remaining-15,420),**kwargs)
    record.append({k:v for k,v in result.items() if k!='samples'})
    (ROOT/'batch_progress.json').write_text(json.dumps({'deadline_utc':'2026-10-19T21:36:00Z','completed':record},indent=2)+'\n')
    if result['reason']=='nonfinite':break

# Use remaining time for additional selected directions rather than infer a
# spectrum from one disturbance. Keep output separate from the earlier study.
if DEADLINE-time.time()>170:
    import azimuthal
    azimuthal.ROOT=ROOT/'long_linear';azimuthal.ROOT.mkdir(exist_ok=True)
    for m in [1,2,3,4]:
        if DEADLINE-time.time()<100:break
        result=azimuthal.run_mode(m,T=96)
        record.append({'type':'linear',**{k:v for k,v in result.items() if k!='samples'}})
        (ROOT/'batch_progress.json').write_text(json.dumps({'deadline_utc':'2026-10-19T21:36:00Z','completed':record},indent=2)+'\n')
(ROOT/'BATCH_FINISHED.json').write_text(json.dumps({'finished_utc':datetime.now(timezone.utc).isoformat(),'items':len(record),'deadline_utc':'2026-10-19T21:36:00Z'},indent=2)+'\n')
print('BATCH FINISHED',flush=True)
