"""Prescribed-geometry diagnostic of the actual linear-coupling internal action.

Known Mathieu/Floquet mechanism; not a new physical effect or joint worldline solution.
"""
from pathlib import Path
import json,time
import numpy as np
from scipy.integrate import solve_ivp
ROOT=Path(__file__).resolve().parent

def floquet(epsilon,Omega=2.,sign=-1,tol=1e-11):
    # h_A=h_S=1, k_A=k_S=2, mu1(A)=A, kappa^2=1+epsilon*cos(Omega*s).
    # Normal coordinates (A +/- S)/sqrt(2) diagonalize the Hessian at all times.
    period=2*np.pi/Omega
    def rhs(t,y):
        M=y.reshape(2,2)
        freq2=2+sign*(1+epsilon*np.cos(Omega*t))
        return (np.array([[0.,1.],[-freq2,0.]])@M).ravel()
    start=time.perf_counter()
    sol=solve_ivp(rhs,(0,period),np.eye(2).ravel(),method='DOP853',rtol=tol,atol=tol/100)
    assert sol.success
    monodromy=sol.y[:,-1].reshape(2,2)
    eig=np.linalg.eigvals(monodromy)
    growth=max(0.,float(np.log(np.max(abs(eig)))/period))
    return dict(epsilon=epsilon,Omega=Omega,sign=sign,tol=tol,period=period,
        minimum_instantaneous_stiffness=1-abs(epsilon),growth_exponent=growth,
        leading_resonant_prediction=epsilon/4 if sign==-1 and Omega==2 else None,
        determinant_error=float(abs(np.linalg.det(monodromy)-1)),
        monodromy=monodromy.tolist(),multipliers=[dict(real=float(z.real),imag=float(z.imag)) for z in eig],
        seconds=time.perf_counter()-start)

if __name__=='__main__':
    rows=[floquet(e,Omega=O,sign=sgn) for e in [.02,.05,.1,.2]
          for O in [1.6,2.,2.4] for sgn in [-1,1]]
    check=floquet(.1,tol=1e-13)
    assert all(r['minimum_instantaneous_stiffness']>0 and r['determinant_error']<1e-9 for r in rows)
    assert all(r['growth_exponent']>0 for r in rows if r['sign']==-1 and r['Omega']==2.)
    (ROOT/'foundations_floquet.json').write_text(json.dumps(dict(runs=rows,refinement=check),indent=2))
    print(json.dumps([r for r in rows if r['sign']==-1 and r['Omega']==2],indent=2))
