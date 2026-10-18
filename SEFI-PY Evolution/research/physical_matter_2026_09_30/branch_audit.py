"""Audited physical stress/current observables, separate from v0.1.

Continuum quadrature and the solver's discrete variational functional are
reported separately; neither is substituted for the other after inspection.
"""
import numpy as np
from scipy.integrate import simpson
from full_z_spectral import operators


def physical(fields,n,h,c,nu,lam):
    u,w,s=fields;p=u+1j*w
    r=np.arange(n+1)*h;z=np.arange(-n,n+1)*h
    integrate=lambda f:float(2*np.pi*simpson(simpson(f,x=z,axis=1)*r,x=r))
    pr,pz=np.gradient(p,h,edge_order=2);sr,sz=np.gradient(s,h,edge_order=2)
    omega=np.sqrt(2);a=c/(2*omega);gamma=np.sqrt(1+a*a);v=a/gamma
    pt=1j*omega*p-a*pz;px=gamma*pz
    st=1j*nu*gamma*s-a*sz;sx=gamma*sz-1j*nu*gamma*v*s
    density=abs(p)**2
    potential=.5*density+.25*density**2+lam/4*s**4+.5*(4*density-3)*s*s
    E=integrate(.5*(abs(pt)**2+abs(st)**2+abs(px)**2+abs(sx)**2+abs(pr)**2+sr**2)+potential-1.75)/gamma
    P=integrate(-np.real(pt.conj()*px+st.conj()*sx))/gamma
    DQ=integrate(np.imag(p.conj()*pt)-omega)/gamma
    Q=integrate(np.imag(s*st))/gamma
    Pred=integrate(np.imag(p.conj()*pz))
    Psurface=integrate(np.imag(pz))
    Kp=.5*integrate(abs(pr)**2+sr**2);Kz=.5*integrate(abs(pz)**2+sz**2)
    U=integrate(.25*(density-1+4*s*s)**2+(lam-16)/4*s**4+(1-nu*nu)/2*s*s)
    G=Kp+Kz+U+c/2*Pred
    H=E-v*P-omega*DQ-nu/gamma*Q
    Rperp=Kz+U+c/2*Pred;Rz=Kp-Kz+U
    lateral=np.pi*(n*h)**2/2*simpson(abs(pr[-1])**2+sr[-1]**2,x=z)
    faces=np.pi*(n*h)*simpson(r*(abs(pz[:,0])**2+sz[:,0]**2+abs(pz[:,-1])**2+sz[:,-1]**2),x=r)
    scale=Kp+Kz+abs(U)+abs(c/2*Pred)
    mass=integrate(s*s)
    peak=np.unravel_index(np.argmax(s),s.shape)
    return {'c':c,'nu':nu,'lambda':lam,'L':n*h,'h':h,'gamma':gamma,'v':v,
        'E_excess':E,'P_physical_z':P,'Delta_Q_phi':DQ,'Q_sigma':Q,
        'P_reduced_bare':Pred,'integral_total_w_Z':Psurface,
        'K_perp':Kp,'K_Z':Kz,'U':U,'G_continuum_quadrature':G,
        'generator_algebra_residual':G-gamma*H,
        'charge_algebra_residual':Q-nu*mass,
        'pohozaev_transverse':Rperp,'pohozaev_longitudinal':Rz,
        'pohozaev_transverse_relative':Rperp/max(scale,1e-30),
        'pohozaev_longitudinal_relative':Rz/max(scale,1e-30),
        'boundary_transverse':float(lateral),'boundary_longitudinal':float(faces),
        'pohozaev_transverse_boundary_corrected':float(Rperp+lateral),
        'pohozaev_longitudinal_boundary_corrected':float(Rz+faces),
        'carrier_norm_squared':mass,'radius_carrier_rms':float(np.sqrt(integrate(r[:,None]**2*s*s)/mass)),
        'radius_carrier_peak':float(r[peak[0]]),'carrier_max':float(s.max()),
        'carrier_axial_rms':float(np.sqrt(integrate(z[None,:]**2*s*s)/mass)),
        'background_parity_defect':max(float(abs(u-u[:,::-1]).max()),float(abs(w+w[:,::-1]).max()),float(abs(s-s[:,::-1]).max()))}


def discrete_functional(fields,n,h,c,nu,lam):
    lap,dz,_=operators(n,h)
    u,w,s=[f[:n,1:-1].ravel() for f in fields];q=u-1
    r=np.arange(n,dtype=float);r[0]=1/8
    weight=2*np.pi*np.repeat(r*h**3,2*n-1)
    # h^3*i = r_i dr dZ; axis weight matches the solver exactly.
    integrate=lambda f:float(np.dot(weight,f))
    K=-.5*integrate(q*(lap@q)+w*(lap@w)+s*(lap@s))
    Pred=integrate(q*(dz@w)-w*(dz@q))
    U=integrate(.25*(u*u+w*w-1+4*s*s)**2+(lam-16)/4*s**4+(1-nu*nu)/2*s*s)
    return {'G_discrete':K+U+c/2*Pred,'K_discrete':K,'U_discrete':U,
            'P_reduced_discrete_variational':Pred,'Q_sigma_discrete':nu*integrate(s*s)}


def diagnostics(fields,n,h,c,nu,lam):
    return {**physical(fields,n,h,c,nu,lam),**discrete_functional(fields,n,h,c,nu,lam)}
