"""Causal conventional estimators. No hidden truth or future data accepted."""
import math
from collections import deque
import numpy as np

R = 1.0 / 32.0

def project(a, b, c):
    """PSD projection of symmetric [[a,b],[b,c]], eigenvalue floor 1e-8."""
    mid=(a+c)*.5; rad=math.hypot((a-c)*.5,b)
    lo=max(1e-8,mid-rad); hi=max(1e-8,mid+rad)
    if rad<1e-15:return lo,0.,lo
    scale=(hi-lo)/(2*rad); center=(hi+lo)*.5
    return center+scale*(a-c)*.5,scale*b,center-scale*(a-c)*.5

def posterior(y, mean, covariance):
    a,b,c=covariance; a+=R; c+=R; det=a*c-b*b
    dx=y[0]-mean[0];dy=y[1]-mean[1]
    ix=(c*dx-b*dy)/det;iy=(a*dy-b*dx)/det
    estimate=(y[0]-R*ix,y[1]-R*iy)
    loglike=-math.log(2*math.pi)-.5*(math.log(det)+dx*ix+dy*iy)
    return estimate,loglike

class Moments:
    def __init__(self, calibration):
        self.mean=np.mean(calibration,axis=0).tolist()
        s=np.cov(calibration,rowvar=False,ddof=0)
        self.obs=(float(s[0,0]),float(s[0,1]),float(s[1,1]))
    def covariance(self):
        a,b,c=self.obs;return project(a-R,b,c-R)
    def update(self,y,rate):
        dx=y[0]-self.mean[0];dy=y[1]-self.mean[1];a,b,c=self.obs
        q=1-rate
        self.obs=(q*(a+rate*dx*dx),q*(b+rate*dx*dy),q*(c+rate*dy*dy))
        self.mean[0]+=rate*dx;self.mean[1]+=rate*dy

class Fixed:
    def __init__(self,calibration):
        self.m=Moments(calibration);self.cov=self.m.covariance()
    def step(self,y):return posterior(y,self.m.mean,self.cov)[0]

class EW:
    def __init__(self,calibration,rate):self.m=Moments(calibration);self.rate=rate
    def step(self,y):
        out=posterior(y,self.m.mean,self.m.covariance())[0]
        self.m.update(y,self.rate);return out

class Rolling:
    def __init__(self,calibration,window):
        self.q=deque(tuple(map(float,y)) for y in calibration[-window:]);self.window=window
        x=np.array(self.q);self.s=x.sum(axis=0);self.ss=x.T@x
    def step(self,y):
        n=len(self.q);mu=self.s/n;s=self.ss/n-np.outer(mu,mu)
        cov=project(float(s[0,0]-R),float(s[0,1]),float(s[1,1]-R))
        out=posterior(y,mu,cov)[0]
        old=np.asarray(self.q.popleft());new=np.asarray(y)
        self.s+=new-old;self.ss+=np.outer(new,new)-np.outer(old,old)
        self.q.append(tuple(new));return out

class Mixture:
    """Standard predictive mixture with fixed-share weight smoothing."""
    def __init__(self,calibration,rates=(.01,.05,.15),share=.005):
        self.experts=[EW(calibration,a) for a in rates];self.weights=np.ones(len(rates))/len(rates);self.share=share
    def step(self,y):
        pairs=[posterior(y,e.m.mean,e.m.covariance()) for e in self.experts]
        ll=np.array([p[1] for p in pairs]);logw=np.log(self.weights)+ll
        w=np.exp(logw-logw.max());w/=w.sum()
        out=w@np.array([p[0] for p in pairs])
        self.weights=(1-self.share)*w+self.share/len(w)
        for e in self.experts:e.m.update(y,e.rate)
        return out

def predict(kind,calibration,observations):
    if kind=='no_recovery':return np.zeros_like(observations)
    if kind=='ml':return np.array(observations,copy=True)
    if kind=='fixed':decoder=Fixed(calibration)
    elif kind.startswith('ew:'):decoder=EW(calibration,float(kind.split(':')[1]))
    elif kind.startswith('rolling:'):decoder=Rolling(calibration,int(kind.split(':')[1]))
    elif kind=='mixture':decoder=Mixture(calibration)
    else:raise ValueError(kind)
    return np.asarray([decoder.step(y) for y in observations])
