"""Matched far-field linear-response benchmark; not a vortex stability run."""
from pathlib import Path
import json, hashlib
import numpy as np
from PIL import Image, ImageDraw, ImageFont

ROOT=Path(__file__).resolve().parent
OMEGA=np.sqrt(2.)

def branches(k,omega):
    A=1+2*omega**2
    # Rationalized lower branch avoids cancellation near k=0.
    disc=np.sqrt(A*A+4*omega**2*k*k)
    low=k*k*(k*k+2)/(k*k+A+disc)
    return np.sqrt(low),np.sqrt(k*k+A+disc)

def matrix(k,omega):
    # State (a,b,a_dot,b_dot), F=1+a+i*b.
    return np.array([[0,0,1,0],[0,0,0,1],[-(k*k+2),0,0,2*omega],[0,-k*k,-2*omega,0]],float)

def integrate(M,y,dt,T):
    n=round(T/dt);h=T/n;y=y.copy()
    for _ in range(n):
        a=M@y;b=M@(y+h*a/2);c=M@(y+h*b/2);d=M@(y+h*c)
        y+=h*(a+2*b+2*c+d)/6
    return y

rows=[];checks=[]
for name,omega in [('static_global_witten_control',0.),('rotating_console_background',OMEGA)]:
    for k in [.05,.1,.25,.5,1.,2.]:
        lo,hi=branches(k,omega); M=matrix(k,omega)
        vals,vecs=np.linalg.eig(M)
        idx=np.argmin(abs(vals-1j*lo));y=vecs[:,idx].astype(complex)
        eigerr=float(abs(vals[idx]-1j*lo))
        T=10.;exact=y*np.exp(1j*lo*T)
        errors=[float(np.max(abs(integrate(M,y,dt,T)-exact))) for dt in [.1,.05]]
        assert eigerr<1e-12 and errors[1]<errors[0]
        rows.append(dict(model=name,k=k,low_frequency=float(lo),high_frequency=float(hi),low_phase_speed=float(lo/k)))
        checks.append(dict(model=name,k=k,eigenvalue_error=eigerr,rk4_dt_01_error=errors[0],rk4_dt_005_error=errors[1]))
assert abs(branches(.0001,OMEGA)[0]/.0001-1/np.sqrt(5))<1e-8
assert np.allclose(branches(np.linspace(0,3,31),0)[0],np.linspace(0,3,31))
result=dict(scope='Dimensionless homogeneous far-field linearization; no nonlinear vortex evolution or empirical validation',omega=float(OMEGA),static_low_k_speed=1.,rotating_low_k_speed=float(1/np.sqrt(5)),speed_reduction_percent=float(100*(1-1/np.sqrt(5))),rows=rows,checks=checks)
(ROOT/'results.json').write_text(json.dumps(result,indent=2))
k=np.linspace(0,3,500)
im=Image.new('RGB',(1400,720),'white'); dr=ImageDraw.Draw(im)
font=ImageFont.truetype('C:/Windows/Fonts/arial.ttf',22)
for panel,ymax,title in [(0,6.,'Co-rotating mode frequencies'),(1,1.1,'Gapless phase speed / light speed')]:
    x0=75+panel*700; y0=540; w=570; h=400
    dr.text((x0,40),title,font=font,fill='black')
    def xy(x,y): return (float(x0+w*x/3),float(y0-h*y/ymax))
    for tick in np.linspace(0,ymax,6):
        yy=xy(0,tick)[1]; dr.line((x0,yy,x0+w,yy),fill='#dddddd'); dr.text((x0-50,yy-10),f'{tick:.2g}',font=font,fill='black')
    for tick in range(4):
        xx=xy(tick,0)[0]; dr.text((xx,y0+10),str(tick),font=font,fill='black')
    dr.text((x0+70,590),'Dimensionless wave number k',font=font,fill='black')
    for om,col in [(0,'#c16a32'),(OMEGA,'#175a7a')]:
        lo,hi=branches(k,om)
        if panel==0:
            dr.line([xy(a,b) for a,b in zip(k,lo)],fill=col,width=4)
            pts=[xy(a,b) for a,b in zip(k,hi)]
            for i in range(0,len(pts)-6,12):dr.line(pts[i:i+6],fill=col,width=3)
        else:dr.line([xy(a,b) for a,b in zip(k[1:],lo[1:]/k[1:])],fill=col,width=4)
dr.text((75,645),'Orange: static global-Witten control | Blue: rotating console background',font=font,fill='black')
dr.text((75,680),'Solid: gapless. Dashed: gapped. Linear far-field test; not nonlinear vortex evolution.',font=font,fill='black')
im.save(ROOT/'dispersion-comparison.png')
print(json.dumps({k:v for k,v in result.items() if k not in ['rows','checks']},indent=2))
print('Verified',len(checks),'mode cases against matrix eigenvalues and RK4 step refinement.')
