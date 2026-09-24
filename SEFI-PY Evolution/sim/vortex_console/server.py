"""Loopback-only research console. No external services or engine imports."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
os.environ.setdefault('OMP_NUM_THREADS','2')
import argparse,json,time,threading,uuid,re,traceback
from pathlib import Path
from http.server import ThreadingHTTPServer,BaseHTTPRequestHandler
from urllib.parse import urlparse,parse_qs
import numpy as np
import physics
import rest_search

ROOT=Path(__file__).resolve().parent;RUNS=ROOT/'runs';RUNS.mkdir(exist_ok=True)
JOBS={};LOCK=threading.Lock();PORT=8765
REFERENCE,BASE=physics.load_profile(ROOT/'data/charged_ring_refined.npz')
BASE_REPORT=json.loads((ROOT/'data/charged_ring_refined.json').read_text())

def dump(path,obj):path.write_text(json.dumps(obj,allow_nan=False),encoding='utf-8')
def safe_id(value):
    if not re.fullmatch('[a-f0-9]{32}',value or ''):raise ValueError('Invalid run identifier')
    return value

class Job:
    model_version=physics.MODEL_VERSION
    def __init__(self,p,action,field=None,values=None,restore=None):
        self.id=uuid.uuid4().hex;self.p=p;self.action=action;self.field=field;self.values=values
        self.dir=RUNS/self.id;self.dir.mkdir();self.state='queued';self.message='Waiting';self.frames=[];self.history=[];self.members=[]
        self.stop=False;self.paused=False;self.want_checkpoint=False;self.ev=None;self.restore=restore;self.error=None;self.result=None
        self.started=time.monotonic();self.finished=None;self.paused_seconds=0.;self.solver_iteration=0;self.residual=None;self.bytes=0;self.current_member=0;self.checkpoint=False
        dump(self.dir/'settings.json',dict(model_version=self.model_version,parameters=p,action=action,sweep_field=field,sweep_values=values,source=p.get('initial','charged_ring_refined')))
    def elapsed(self):return (self.finished or time.monotonic())-self.started-self.paused_seconds
    def guard(self):
        if self.stop:raise InterruptedError('Stopped by user')
        if self.elapsed()>self.p['wall_seconds']:raise InterruptedError('Runtime limit reached; checkpoint saved')
    def cancelled(self):self.guard();return False
    def status(self):
        d=self.history[-1] if self.history else {}
        sim=d.get('time',0);elapsed=self.elapsed();fraction=sim/self.p.get('duration',1)
        return dict(id=self.id,model_version=self.model_version,state=self.state,message=self.message,parameters=self.p,action=self.action,
          frames=len(self.frames),history=self.history[-201:],members=self.members,member=self.current_member,
          elapsed=round(elapsed,2),eta=round(elapsed*(1-fraction)/fraction,1) if fraction>0 and self.state=='running' else None,
          residual=self.residual,iteration=self.solver_iteration,result=self.result,error=self.error,checkpoint=self.checkpoint,
          storage_bytes=self.bytes)
    def save_frame(self,f):
        f['model_version']=self.model_version;f['member']=self.current_member;f['parameters']=self.p.copy();data=json.dumps(f,allow_nan=False).encode()
        if self.bytes+len(data)+1000000>self.p['storage_mb']*1024**2:raise InterruptedError('Storage limit reached; checkpoint saved')
        name=f'frame_{len(self.frames):04d}.json';(self.dir/name).write_bytes(data);self.bytes+=len(data);self.frames.append(name)
    def save_checkpoint(self):
        if self.ev is None:return
        e=self.ev
        temp=self.dir/'checkpoint.tmp.npz'
        np.savez_compressed(temp,model_version=physics.MODEL_VERSION,y=e.y,time=e.time,steps=e.steps,q0=e.q0,parameters=json.dumps(self.p),
             **{k:self.profile[k] for k in ['r','z','u','w','s']})
        temp.replace(self.dir/'checkpoint.npz');self.checkpoint=True;self.want_checkpoint=False
    def tick(self,i,res):self.solver_iteration=i;self.residual=res;self.message=f'Profile iteration {i}; residual {res:.3g}'
    def get_profile(self,p):
        same=all(abs(p[k]-BASE[k])<1e-12 for k in ['c','nu','lam','N','L','n']) and abs(p['seed']-5.36)<.02
        if same and self.action!='solve':
            self.residual=BASE_REPORT['residual'];return {k:v.copy() for k,v in REFERENCE.items()}
        self.state='solving';self.message='Solving both fields';f,result=physics.solve(p,REFERENCE,self.tick,self.cancelled)
        self.result=result;self.residual=result['residual']
        if not result['converged'] or not result['has_ring'] or result['carrier_peak']<1e-5:
            pp,ss,z=physics.full_fields(f);self.save_frame(physics.frame(f['r'],z,pp,ss,dict(time=0,residual=self.residual)))
            raise ValueError('No converged charged ring: '+result['reason']+'. Trial is shown; evolution was not started.')
        return f
    def run(self):
        try:
            values=self.values if self.values is not None else [None]
            for member,value in enumerate(values):
                self.current_member=member
                if value is not None:self.p[self.field]=value
                if self.restore:
                    d=np.load(self.restore,allow_pickle=False);self.profile={k:d[k] for k in ['r','z','u','w','s']}
                else:self.profile=self.get_profile(self.p)
                np.savez_compressed(self.dir/f'profile_{member:02d}.npz',model_version=physics.MODEL_VERSION,parameters=json.dumps(self.p),**self.profile)
                if self.action=='solve':
                    pp,ss,z=physics.full_fields(self.profile);rs=physics.roots(self.profile)
                    self.save_frame(physics.frame(self.profile['r'],z,pp,ss,dict(time=0,radius=rs[0],residual=self.residual)))
                    self.state='completed';self.message='Converged numerical profile; stability untested';return
                self.ev=physics.Evolution(self.profile,self.p)
                if self.restore:
                    self.ev.y=d['y'].copy();self.ev.time=float(d['time']);self.ev.steps=int(d['steps']);self.ev.q0=float(d['q0'])
                self.state='running';self.message='Integrating the coupled field equations';last_save=-1e9
                while True:
                    self.guard()
                    if self.paused:
                        self.state='paused';t=time.monotonic()
                        while self.paused and not self.stop:
                            if self.want_checkpoint:self.save_checkpoint()
                            time.sleep(.03)
                        self.paused_seconds+=time.monotonic()-t;self.state='running';self.guard()
                    if self.want_checkpoint:self.save_checkpoint()
                    e=self.ev
                    if e.time-last_save>=self.p['save_every']-1e-9 or e.time>=self.p['duration']-1e-9:
                        diag=e.diagnostics();diag['member']=member;self.history.append(diag)
                        self.save_frame(physics.frame(e.r,e.z,e.y[0],e.y[1],diag));last_save=e.time
                    if e.time>=self.p['duration']-1e-9:break
                    e.step(min(self.p['dt'],self.p['duration']-e.time))
                self.save_checkpoint();self.members.append(dict(index=member,parameters=self.p.copy(),final=self.history[-1]))
            self.state='completed';self.message='Computed evolution complete. Limited axisymmetric test, not a stability proof.'
        except InterruptedError as exc:self.state='stopped';self.message=str(exc)
        except Exception as exc:self.state='failed';self.error=str(exc);self.message=str(exc)
        finally:
            self.finished=time.monotonic()
            try:self.save_checkpoint();dump(self.dir/'result.json',self.status())
            except Exception as exc:self.error=(self.error or '')+'; could not save: '+str(exc)

def rest_frame(search,iteration):
    d=search.diagnostics();d['iteration']=iteration
    return {**physics.frame(search.r,search.z,search.y[0]+1j*search.y[1],search.y[2],d),'kind':'relaxation','source':'Fixed-charge relaxation; not physical time evolution'}

class RestJob(Job):
    model_version=rest_search.MODEL_VERSION
    def save_checkpoint(self):
        if not hasattr(self,'search'):return
        temp=self.dir/'rest_checkpoint.tmp.npz'
        np.savez_compressed(temp,model_version=self.model_version,parameters=json.dumps(self.p),y=self.search.y,r=self.search.r,z=self.search.z,iteration=self.solver_iteration)
        temp.replace(self.dir/'rest_checkpoint.npz');self.checkpoint=True;self.want_checkpoint=False
    def run(self):
        try:
            p=self.p;self.search=rest_search.Search(**{k:p[k] for k in ['L','h','Q','lam','N']})
            if self.restore:
                with np.load(self.restore,allow_pickle=False) as d:self.search.y=d['y'].copy();self.solver_iteration=int(d['iteration'])
            elif p['initial']!='ring':rest_search.load_initial(self.search,ROOT/f'data/fixed_charge_{p["initial"]}.npz')
            self.state='running';self.message='Fixed-charge relaxation — iteration count is not physical time'
            while True:
                self.guard()
                if self.paused:
                    self.state='paused';t=time.monotonic()
                    while self.paused and not self.stop:
                        if self.want_checkpoint:self.save_checkpoint()
                        time.sleep(.03)
                    self.paused_seconds+=time.monotonic()-t;self.state='running';self.guard()
                if self.want_checkpoint:self.save_checkpoint()
                i=self.solver_iteration;res=float(abs(self.search.force(self.search.y)).max());self.residual=res
                done=res<1e-5 or i>=p['steps']
                if not self.history or i%p['every']==0 or done:
                    f=rest_frame(self.search,i);self.history.append(f['diagnostics']);self.save_frame(f)
                if done:
                    self.state='completed';self.result=self.search.diagnostics()
                    self.message=('Residual tolerance reached; candidate only. ' if res<1e-5 else 'Iteration limit reached; unconverged trial. ')+('No original-field midplane vortex core. ' if not self.result['core_radii'] else 'Midplane core crossings present; topology needs verification. ')+'Stability and particle mass are unproved.'
                    break
                self.search.descent_step();self.solver_iteration+=1
        except InterruptedError as exc:self.state='stopped';self.message=str(exc)
        except Exception as exc:self.state='failed';self.error=str(exc);self.message=str(exc)
        finally:
            self.finished=time.monotonic()
            try:self.save_checkpoint();dump(self.dir/'result.json',self.status())
            except Exception as exc:self.error=(self.error or '')+'; save failed: '+str(exc)

def new_rest_job(body):
    p=rest_search.validate(body.get('parameters',{}));restore=None
    if body['action']=='restore_rest':
        restore=RUNS/safe_id(body.get('checkpoint'))/'rest_checkpoint.npz'
        with np.load(restore,allow_pickle=False) as d:
            if str(d.get('model_version',''))!=rest_search.MODEL_VERSION:raise ValueError('Incompatible relaxation checkpoint')
            old=json.loads(str(d['parameters']));old.update({k:p[k] for k in ['steps','wall_seconds','storage_mb']});p=rest_search.validate(old)
            if p['steps']<=int(d['iteration']):raise ValueError('Maximum iteration must exceed checkpoint iteration')
    with LOCK:
        if any(j.state in ['queued','solving','running','paused'] for j in JOBS.values()):raise ValueError('Stop or finish the active calculation first')
        if sum(f.stat().st_size for f in RUNS.rglob('*') if f.is_file())>500*1024**2:raise ValueError('Archive local runs before starting more: 500 MB limit reached')
        j=RestJob(p,body['action'],restore=restore);JOBS[j.id]=j
    threading.Thread(target=j.run,daemon=True).start();return j

def new_job(body):
    action=body.get('action','evolve')
    if action in ['relax','restore_rest']:return new_rest_job(body)
    if action not in ['solve','evolve','sweep','restore']:raise ValueError('Unknown action')
    p=physics.validate(body.get('parameters',{}),action!='solve');field=None;values=None;restore=None
    if action=='sweep':
        field=body.get('field');values=body.get('values')
        if field not in ['epsilon','c','nu'] or not isinstance(values,list) or not 1<=len(values)<=5:raise ValueError('Sweep requires epsilon, c or nu and one to five values')
        for value in values:
            candidate=physics.validate({**p,field:value},True)
        values=[float(v) for v in values]
    if action=='restore':
        ident=safe_id(body.get('checkpoint'));restore=RUNS/ident/'checkpoint.npz'
        if not restore.exists():raise ValueError('Checkpoint not found')
        with np.load(restore,allow_pickle=False) as d:
            if str(d.get('model_version',''))!=physics.MODEL_VERSION:raise ValueError('Checkpoint model version differs from this console; explicit migration is required')
            original=json.loads(str(d['parameters']));original['duration']=p['duration'];original['wall_seconds']=p['wall_seconds'];original['storage_mb']=p['storage_mb']
            p=physics.validate(original,True)
            if p['duration']<=float(d['time']):raise ValueError('End time must be later than the checkpoint time')
    with LOCK:
        if any(j.state in ['queued','solving','running','paused'] for j in JOBS.values()):raise ValueError('Stop or finish the active calculation first')
        if sum(f.stat().st_size for f in RUNS.rglob('*') if f.is_file())>500*1024**2:raise ValueError('Saved run storage exceeds 500 MB. Archive runs before starting more.')
        j=Job(p,action,field,values,restore);JOBS[j.id]=j
    threading.Thread(target=j.run,daemon=True).start();return j

class Handler(BaseHTTPRequestHandler):
    def log_message(self,*args):pass
    def respond(self,obj,status=200):
        data=json.dumps(obj,allow_nan=False).encode();self.send_response(status);self.send_header('Content-Type','application/json');self.send_header('Content-Length',str(len(data)));self.send_header('Cache-Control','no-store');self.end_headers();self.wfile.write(data)
    def check_host(self):
        return self.headers.get('Host') in [f'127.0.0.1:{PORT}',f'localhost:{PORT}']
    def do_GET(self):
        if not self.check_host():self.respond({'error':'Invalid host'},403);return
        u=urlparse(self.path);q=parse_qs(u.query)
        try:
            if u.path=='/api/rest_baseline':
                name=q.get('name',['Q1000_wide'])[0]
                if name not in rest_search.SAVED:raise ValueError('Unknown rest reference')
                with np.load(ROOT/f'data/fixed_charge_{name}.npz',allow_pickle=False) as d:
                    p=json.loads(str(d['parameters']));s=rest_search.Search(**p);s.y=d['y'].copy()
                report=json.loads((ROOT/f'data/fixed_charge_{name}.json').read_text());f=rest_frame(s,report['iterations']);f.update(parameters={**rest_search.DEFAULTS,**p},model_version=rest_search.MODEL_VERSION,report=report)
                self.respond(f);return
            if u.path=='/api/baseline':
                name=q.get('name',['refined'])[0]
                if name not in ['refined','trial','wide']:raise ValueError('Unknown profile')
                f,p=physics.load_profile(ROOT/f'data/charged_ring_{name}.npz');report=json.loads((ROOT/f'data/charged_ring_{name}.json').read_text())
                pp,ss,z=physics.full_fields(f);data=physics.frame(f['r'],z,pp,ss,dict(time=0,charge=report['carrier_charge'],winding=1,residual=report['residual'],radius=report['core_radii'][0]));data.update(parameters=p,source='Saved numerical profile',report=report)
                self.respond(data);return
            if u.path=='/api/jobs':
                checkpoints=[p.parent.name for p in RUNS.glob('*/checkpoint.npz')]
                self.respond(dict(jobs=[j.status() for j in JOBS.values()],checkpoints=checkpoints,rest_checkpoints=[p.parent.name for p in RUNS.glob('*/rest_checkpoint.npz')]));return
            if u.path in ['/api/status','/api/frame','/api/checkpoint','/api/rest_checkpoint','/api/settings']:
                ident=safe_id(q.get('id',[''])[0]);j=JOBS.get(ident)
                if u.path=='/api/status':
                    if not j:raise ValueError('Run not in this server session')
                    self.respond(j.status());return
                if u.path=='/api/frame':
                    index=int(q.get('index',['0'])[0]);files=j.frames if j else [x.name for x in sorted((RUNS/ident).glob('frame_*.json'))]
                    if not 0<=index<len(files):raise ValueError('Frame index outside saved range')
                    self.respond(json.loads((RUNS/ident/files[index]).read_text()));return
                path=RUNS/ident/('rest_checkpoint.npz' if u.path=='/api/rest_checkpoint' else 'checkpoint.npz' if u.path.endswith('checkpoint') else 'settings.json')
                self.file(path,'application/octet-stream',download=path.name);return
            if u.path=='/research-note':self.file(ROOT/'FIXED_CHARGE_RESEARCH.md','text/plain; charset=utf-8');return
            files={'/':'index.html','/app.js':'app.js','/rest.js':'rest.js','/style.css':'style.css'}
            if u.path not in files:self.respond({'error':'Not found'},404);return
            name=files[u.path];self.file(ROOT/'static'/name,{'html':'text/html','js':'text/javascript','css':'text/css'}[name.split('.')[-1]])
        except (ValueError,KeyError,FileNotFoundError) as e:self.respond({'error':str(e)},400)
    def file(self,path,mime,download=None):
        data=path.read_bytes();self.send_response(200);self.send_header('Content-Type',mime);self.send_header('Content-Length',str(len(data)))
        self.send_header('X-Content-Type-Options','nosniff')
        if download:self.send_header('Content-Disposition',f'attachment; filename="{download}"')
        self.end_headers();self.wfile.write(data)
    def do_POST(self):
        origin=self.headers.get('Origin')
        if not self.check_host() or self.headers.get('X-Vortex-Console')!='1' or (origin and origin not in [f'http://127.0.0.1:{PORT}',f'http://localhost:{PORT}']):
            self.respond({'error':'Only local console requests are accepted'},403);return
        try:
            length=int(self.headers.get('Content-Length','0'))
            if not 0<length<=20000:raise ValueError('Invalid request size')
            body=json.loads(self.rfile.read(length))
            if self.path=='/api/start':j=new_job(body);self.respond(j.status());return
            if self.path=='/api/control':
                ident=safe_id(body.get('id'));j=JOBS[ident];cmd=body.get('command')
                if cmd=='pause':
                    if j.state!='running':raise ValueError('Pause is available during evolution')
                    j.paused=True
                elif cmd=='resume':j.paused=False
                elif cmd=='stop':j.stop=True;j.paused=False
                elif cmd=='checkpoint':
                    if j.ev is None and not hasattr(j,'search'):raise ValueError('Calculation has not started')
                    if j.state in ['completed','stopped','failed']:j.save_checkpoint()
                    else:j.want_checkpoint=True
                else:raise ValueError('Unknown command')
                self.respond(j.status());return
            self.respond({'error':'Not found'},404)
        except (ValueError,KeyError,TypeError,FileNotFoundError) as e:self.respond({'error':str(e)},400)

def main():
    global PORT
    parser=argparse.ArgumentParser();parser.add_argument('--port',type=int,default=8765);a=parser.parse_args();PORT=a.port
    server=ThreadingHTTPServer(('127.0.0.1',PORT),Handler)
    print(f'Vortex Research Console: http://127.0.0.1:{PORT}',flush=True)
    try:server.serve_forever()
    except KeyboardInterrupt:pass
    finally:
        for j in JOBS.values():j.stop=True;j.paused=False
        server.server_close()
if __name__=='__main__':main()

