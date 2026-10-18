# Linearization and spectral-stability formulation

**Date:** 2026-09-30  
**Status:** exact linearization checkpoint; spectrum not yet computed  
**Depends on:** \`action_traveling_audit.md\`  
**Baseline implementation preserved:** \`coupled-ring-axisymmetric-v0.1\`

## Scope

This note derives the second variation and the full time-dependent linearized equations about a stationary traveling profile. It deliberately distinguishes a constrained Hessian test from a dynamical spectral-stability calculation.

## 1. Reduced constrained Hessian

Write the stationary profile as
\[
\psi=u+iw,\qquad \sigma=a+ib,
\]
with computed carrier profile \(a=s,b=0\). The reduced functional is
\[
\mathcal G_c=\int d^3X\left[\frac12(|\nabla\psi|^2+|\nabla\sigma|^2)+\frac{N^2}{2r^2}|\sigma|^2+U_\nu+\frac c2\operatorname{Im}(\psi^*\partial_Z\psi)\right].
\]
At \(a=s,b=0\), define
\[
A=u^2+w^2-1+4s^2,\qquad B=4(u^2+w^2)-3-\nu^2+\lambda s^2.
\]
For \(q=(\delta u,\delta w,\delta a,\delta b)^T\), the constrained Hessian is
\[
\mathcal H=
\begin{pmatrix}
-\Delta+A+2u^2 & 2uw+c\partial_Z & 8us & 0\\
2uw-c\partial_Z & -\Delta+A+2w^2 & 8ws & 0\\
8us & 8ws & -\Delta+N^2/r^2+B+2\lambda s^2 & 0\\
0&0&0&-\Delta+N^2/r^2+B
\end{pmatrix}.
\]
Under localized/infinite-domain boundary conditions, \(\partial_Z^\dagger=-\partial_Z\), so the off-diagonal traveling terms are mutual adjoints.

The \(\delta b\) channel is absent from the current stationary Newton solve because the carrier is gauge-fixed real. Therefore the Newton Jacobian is a restricted three-real-component object, not a full dynamical stability operator.

For the \((u,w,s)\) sector, the negative Hessian has the same local coefficients as the current \`physics.py\` Newton Jacobian: \(1-3u^2-w^2-4s^2\), \(1-u^2-3w^2-4s^2\), \(3+\nu^2-4(u^2+w^2)-3\lambda s^2\), and cross couplings \(-2uw,-8us,-8ws\), with the traveling derivative supplied by the ring operator. This is a structural action-to-code check, not a numerical spectral result.

## 2. Exact time-dependent envelope equations

Use
\[
T=t,\quad Z=\gamma(z-vt),\quad
\Phi=e^{i\Omega t}p(T,r,Z),\quad
\Sigma=e^{i\nu\gamma(t-vz)}s(T,r,Z),
\]
with \(\Omega=\sqrt2\), \(c=2\Omega\gamma v\). Direct substitution gives
\[
p_{TT}=2\gamma v\,p_{TZ}-2i\Omega p_T+\Delta p+icp_Z+(1-|p|^2-4|s|^2)p,
\]
\[
s_{TT}=2\gamma v\,s_{TZ}-2i\nu\gamma s_T+\Delta s-\frac{N^2}{r^2}s+(3+\nu^2-4|p|^2-\lambda|s|^2)s.
\]
For \(N=0,\lambda=25\), these coefficients reproduce the current \`physics.py\` evolution equations.

## 3. Full linearization

Let \(p=p_0+\xi\), \(s=s_0+\eta\), with \(p_0=u+iw\) and \(s_0=s\) real. Then
\[
\delta[(1-|p|^2-4|s|^2)p]
=(1-2|p_0|^2-4s^2)\xi-p_0^2\xi^*-4p_0s(\eta+\eta^*),
\]
\[
\delta[(3+\nu^2-4|p|^2-\lambda|s|^2)s]
=(3+\nu^2-4|p_0|^2-2\lambda s^2)\eta-\lambda s^2\eta^*-4s(p_0^*\xi+p_0\xi^*).
\]
Hence
\[
\xi_{TT}=2\gamma v\,\xi_{TZ}-2i\Omega\xi_T+\Delta\xi+ic\xi_Z+(1-2|p_0|^2-4s^2)\xi-p_0^2\xi^*-4p_0s(\eta+\eta^*),
\]
\[
\eta_{TT}=2\gamma v\,\eta_{TZ}-2i\nu\gamma\eta_T+\Delta\eta-\frac{N^2}{r^2}\eta+(3+\nu^2-4|p_0|^2-2\lambda s^2)\eta-\lambda s^2\eta^*-4s(p_0^*\xi+p_0\xi^*).
\]

## 4. Quadratic spectral problem

For normal modes proportional to \(e^{\sigma T}\), conjugate coupling requires independent Bogoliubov components, e.g.
\[
Y=(\xi,\bar\xi,\eta,\bar\eta)^T.
\]
The spectral problem has quadratic form
\[
[\sigma^2I+\sigma\mathcal C+\mathcal K]Y=0,
\]
where \(\mathcal C\) contains mixed-advection and gyroscopic terms and \(\mathcal K\) contains spatial and Bogoliubov/Hessian couplings. It can be converted to a first-order generalized eigenproblem by introducing \(V=\sigma Y\).

A mode with \(\operatorname{Re}\sigma>0\) is a linear spectral instability. Absence of such modes in a finite discretization is only numerical evidence and requires grid/domain refinement.

## 5. Nonaxisymmetric perturbations

A full perturbation can be decomposed in azimuthal sectors \(e^{im\varphi}\). Each sector introduces the appropriate cylindrical centrifugal term. For carrier background winding \(N\), Bogoliubov components generally carry shifted angular indices \(N\pm m\); these terms must be derived consistently before implementation. Baseline \(N=0\) is the simplest first spectral calculation.

Previous targeted \(m=1,m=2\) time-domain probes do not replace an eigenvalue search because they sample particular initial data and can miss unstable eigenvectors or slow growth.

## 6. Falsification criterion

The traveling branch fails the stability stage if a reproducible eigenvalue with \(\operatorname{Re}\sigma>0\) persists under appropriate grid/domain refinement in any physically allowed perturbation sector.

Failure to find such an eigenvalue is not a proof of nonlinear stability; it is only spectral evidence over the searched sectors and resolutions.

## 7. Claims not made

No spectrum has been computed in this checkpoint. This derivation does not establish stability, matter, particle identity, spin, electric charge, quantization, or correspondence with an observed particle.
