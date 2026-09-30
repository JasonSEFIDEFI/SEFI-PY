# Three-track identities and characteristic audit

This uses the canonical action in the September 30 checkpoint, with signature
\((+---)\), \(\Omega=\sqrt2,\lambda=25,N=0\). Model assumptions are retained;
matter, gravity, and infinite-domain existence are not assumed.

## Audited laboratory quantities

Set \(p=u+iw\), \(a=c/(2\Omega)=\gamma v\),
\(\gamma=\sqrt{1+a^2}\), \(v=a/\gamma\). At laboratory time zero,

\[
\Phi_t=i\Omega p-a p_Z,\quad \Phi_z=\gamma p_Z,\quad
\Sigma_t=i\nu\gamma s-a s_Z,\quad
\Sigma_z=\gamma s_Z-i\nu\gamma v s.
\]

Laboratory volume is \(d^3x=d^3X/\gamma\). Directly integrate the canonical
stress and phase currents:

\[
E_{\rm exc}=\int d^3x\left[
\tfrac12(|\Phi_t|^2+|\Sigma_t|^2+|\nabla\Phi|^2+|\nabla\Sigma|^2)
+V-\tfrac74\right],
\]
\[
P^z=-\int d^3x\,{\rm Re}(\Phi_t^*\Phi_z+\Sigma_t^*\Sigma_z),\quad
\Delta Q_\Phi=\int d^3x[{\rm Im}(\Phi^*\Phi_t)-\Omega],\quad
Q_\Sigma=\nu\int d^3X\,s^2.
\]

These quantities are not GP energy or GP impulse. In particular,
\[
P^z=a\int d^3X(|p_Z|^2+s_Z^2+\nu^2s^2)-\Omega P_{\rm red},\quad
P_{\rm red}=\int d^3X\,{\rm Im}(p^*p_Z).
\]
Pointwise algebra gives
\[
G=\gamma[E_{\rm exc}-vP^z-\Omega\Delta Q_\Phi-(\nu/\gamma)Q_\Sigma]
=K_\perp+K_Z+U+(c/2)P_{\rm red}.
\]
The code checks this equality independently of profile convergence. Negative
background-subtracted energy is retained as measured; it is not an inferred
negative particle mass or a positivity theorem for the relative generator.

## Boundary-aware Pohozaev and first-law checks

For an exact finite cylinder with \(p=1,s=0\) on its outer surface, the
reduced canonical stress gives
\[
R_\perp=K_Z+U+(c/2)P_{\rm red}=-B_\perp,\qquad
R_Z=K_\perp-K_Z+U=-B_Z,
\]
\[
B_\perp=\frac{\pi L^2}{2}\int_{-L}^{L}(|p_r(L,Z)|^2+s_r(L,Z)^2)dZ,
\]
\[
B_Z=\pi L\int_0^L r\left[
|p_Z(r,L)|^2+s_Z(r,L)^2+
|p_Z(r,-L)|^2+s_Z(r,-L)^2\right]dr.
\]
To derive this, use \(T_{ij}=\partial L_{\rm red}/\partial q_i\ q_j-
\delta_{ij}L_{\rm red}\), \(\partial_i T_{ij}=0\), with real \(q=(u,w,s)\).
Its transverse trace is \(-2R_\perp\), and its longitudinal trace is \(-R_Z\).
Tangential derivatives vanish on each fixed Dirichlet face; the advective
term cancels in the longitudinal trace and normal stress. Corners have zero
surface measure. The uncorrected identities must tend to zero as the box
is removed; boundary-corrected identities separately measure discretization.

At exact critical points on a fixed comoving box the envelope identities are
\[
\partial_cG=\tfrac12P_{\rm red},\qquad \partial_\nu G=-Q_\Sigma.
\]
Write \(\mu=\nu/\gamma\), \(H=G/\gamma\), and differentiate
\(G=\gamma(E-vP-\Omega\Delta Q_\Phi-\mu Q_\Sigma)\). Since
\(dv/dc=(2\Omega\gamma^3)^{-1}\), \(d\gamma/dc=v/(2\Omega)\), the
finite-box first-law defect for a c tangent is
\[
F_c=\partial_cE-v\partial_cP-\Omega\partial_c\Delta Q_\Phi
-(\nu/\gamma)\partial_cQ_\Sigma
=-\frac{v}{2\Omega\gamma^2}R_Z.
\]
For a nu tangent at fixed c, \(F_\nu=0\). Both defects vanish in the localized
continuum/domain limit if the relevant limits exist. Charge variations cannot
be dropped: \(dE/dP=v\) is not the general branch identity.

The solver is exactly variational for a distinct graph functional:
weights \(2\pi i h^3\), axis \(i=1/8\), homogeneous perturbations \(q=u-1,w,s\),
gradient energy \(-\frac12\langle q,\Delta_hq\rangle\), and
\[
P_{{\rm red},h}=\langle (u-1)D_Zw-wD_Z(u-1)\rangle.
\]
Its variation matches the committed finite-difference equations. The
subtracted total derivative is zero in the continuum with fixed outer w=0,
but is not exactly zero for a bare interior-node derivative sum. Its discrete
envelope tests are therefore reported separately from continuum Simpson
quadrature. Physical E, P, Q and raw Pohozaev tests are never redefined to
force agreement.

## Full coupled microscopic principal symbol

In Cartesian real fields \(F=(f_1,f_2,g_1,g_2)\), the linearized equations are
\[
[\eta^{\mu\nu}\partial_\mu\partial_\nu I_4+M(F_0)]\delta F=0,
\]
\[
M=
\begin{pmatrix}
(1+A+4B)I_2+2ff^T & 8fg^T\\
8gf^T & (\lambda B+4A-3)I_2+2\lambda gg^T
\end{pmatrix},\quad A=f^Tf,\ B=g^Tg.
\]
All coupling is lower differential order. Thus
\(\mathcal P(k)=(\eta^{\mu\nu}k_\mu k_\nu)I_4\) and
\(\det\mathcal P=(\eta^{\mu\nu}k_\mu k_\nu)^4\).
The exact characteristics are the already-assumed Minkowski null cone at
every background, including cores where polar variables are singular.
Rotating frames add gyroscopic first derivatives, not a new principal cone.
In \(T=t,Z=\gamma(z-vt)\), the same symbol is
\[
[\tau^2-2a\tau k_Z-k_Z^2-|k_\perp|^2]I_4,
\]
which is simply the coordinate transform of Minkowski; physical front speed
is one.

## Infrared cones are sector-dependent

At \(A=1,B=0,\Omega^2=2\), the exact Phi branches obey
\[
\omega_\pm^2=k^2+5\pm\sqrt{25+8k^2}.
\]
The gapless branch has \(c_s^2=1/5\) at small k, while both branches tend
to front speed one at high k. The Sigma field has laboratory dispersion
\(\omega_\Sigma^2=k^2+1\); its phase is undefined on the zero-condensate
background. It is not another massless field in the Phi phase metric.

For a homogeneous two-condensate branch, neglecting amplitude gradients gives
\[
X=(\partial\theta)^2=1+A+4B,\quad
Y=(\partial\chi)^2=-3+4A+\lambda B,\quad D=\lambda-16,
\]
\[
P(X,Y)=\frac{\lambda(X-1)^2-8(X-1)(Y+3)+(Y+3)^2}{4D}.
\]
Elimination requires \(A,B>0\), positive susceptibility, and a derivative
expansion below amplitude gaps. With \(w^\mu=\partial^\mu\theta\) and
\(z^\mu=\partial^\mu\chi\), the coupled phase symbol is
\[
\begin{pmatrix}
A k^2+(2\lambda/D)(w\cdot k)^2 &
-(8/D)(w\cdot k)(z\cdot k)\\
-(8/D)(w\cdot k)(z\cdot k) &
B k^2+(2/D)(z\cdot k)^2
\end{pmatrix}.
\]
In a common rest frame its time matrix is
\[
T={\rm diag}(A,B)+\frac2D
\begin{pmatrix}\lambda\Omega^2&-4\Omega\nu\\-4\Omega\nu&\nu^2\end{pmatrix},
\]
and its spatial matrix is \({\rm diag}(A,B)\). The two sound speeds are the
generalized eigenvalues of these matrices. For A=1,B=0.1, lambda=25 they are
0.04959603434 and 0.33918901239 (squared), reproduced by the exact coupled
dispersion at small k. This illustrative background is not the localized
ring's asymptotic state. On that state, substituting X=2,Y=0.09 into the
two-condensate algebra gives B<0, explicitly invalidating that branch of the
elimination. Strong gradients and zero amplitudes also invalidate extending
the phase-only metric through the ring core.

A generic shared infrared phase metric fails this comparison. A common
microscopic Minkowski cone is present by assumption, not dynamically emerged.
No gravitational field equation, universal coupling, or matter identification
has been derived. This is a negative result for interpreting this model's
phase metric as universal geometry, not a theorem about other models.
