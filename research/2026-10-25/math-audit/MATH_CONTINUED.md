# Mathematical audit and continuation

25 October 2026. SEFI research context: Jason Duran Dutton.

This audits the principal claims developed in this chat and the cited foundational action. It does not certify every manuscript or equation anywhere in the repositories. No experimental validation or priority claim is made.

## Audit map

| Strand | Checked result | Limit |
|---|---|---|
| Foundational action | Internal Euler–Lagrange signs, autonomous Origin energy, square completion and general Schur complement | Fixed-geometry internal stability is not coupled-worldline stability |
| Projection | Orthogonal projection and Gram correction agree with conventional constrained estimation | Requires the stated metric/frame; conserved oscillator energy does not guarantee phase accuracy |
| QEC | Matched-prior geometric scores equal conventional likelihood; smoothing satisfies conventional quadratic stationarity | No unique advantage follows from relabeling likelihood or regularization |
| Nuclear regression | Units, disjoint data partitions and conventional GP algebra pass selected tests | Regression does not establish new nuclei; heavy-proton extrapolation remains poor |
| Voxel | Plane-wave cancellation and voxel accumulation pass selected tests | Not a calibrated propagation or multipoint trapping solution |
| Two-field model | Binding threshold, induced metric and leading scaling survive | Charge fold, binding and dynamical stability differ |

Sources: output/theory-foundations/foundations.tex; SEFI-PY/research/2026-10-19/vortex-extended/vacuum_spherical.py; the Oct 23 critical-binding notebook; selected QEC, geometric, nuclear and Voxel test logs saved here.

Corrections and scope:

1. Boundary-value convergence alone does not identify a physical branch. Check frequency window, nontriviality, monotonicity, domain and virial identities.
2. The initial L=320 low-amplitude run failed a stronger virial check. Direct initialization at L=480 converged to negative alpha. Both logs remain; continuation from an accepted shorter-domain profile repaired initialization.
3. The old extra Sigma-phase block belongs to a complex-Sigma extension. This model has real Sigma: two real amplitudes plus the charged field's phase are the appropriate directions.
4. The nonlocal fixed-charge second variation includes a positive rank-one term. Dropping it changes radial stability.
5. The number 16 comes from the chosen dimensionless potential coefficients, not an identified constant of nature.

## Exact general quartic criterion

For \(a,c,\lambda,v^2>0\) and \(b\ge0\), with \(m^2=b+cv^2>0\), consider
\[
 U=\frac a4f^4+\frac b2f^2+\frac c2f^2s^2+
 \frac\lambda4(s^2-v^2)^2.
\]
The repository has \(a=b=1,c=4,v^2=3/\lambda\). Exactly,
\[
 U-\frac12m^2f^2=
 \frac\lambda4\left(s^2-v^2+\frac c\lambda f^2\right)^2+
 \frac{a\lambda-c^2}{4\lambda}f^4.
\]
With \(I=\int f^2\), positive gradients, and
\(E=T+\int U+Q^2/(2I)\), if \(a\lambda\ge c^2\) then
\[
 E\ge\frac12m^2I+\frac{Q^2}{2I}\ge m|Q|.
\]
For a regular localized stationary state with \(|\omega|\le m\), \(U_\omega=U-\omega^2f^2/2\ge0\).
The three-dimensional dilation identity \(T+3\int U_\omega=0\) then forces triviality, subject to the stated regularity/decay assumptions.

When \(a\lambda<c^2\), constrained minimization over \(s^2\ge0,f^2>0\) instead gives
\[
 \omega_{\min}^2=b+v^2\sqrt{a\lambda}<m^2,\qquad
 f_*^2=v^2\sqrt{\lambda/a},\quad s_*=0.
\]
Consequently \(E\ge\omega_{\min}|Q|\), and large-core trial states approach this ratio from above. This is not a complete existence/compactness theorem. The natural control ratio is \(c^2/(a\lambda)\), explaining the special value \(\lambda_c=16\).

## Two different thresholds

Set \(\delta=16-\lambda,\ \omega^2=m_\lambda^2-\delta\alpha,\ \rho=\sqrt\delta R,\ m_c=\sqrt{7/4}\).
The critical valley and reduced potential are
\[
 f_0=\frac{\sqrt3}{2}\sin\theta,\quad s_0=\frac{\sqrt3}{4}\cos\theta,
 \quad G=\frac3{16}(1+3\cos^2\theta),
\]
\[
 V_\alpha=\frac{3\alpha}{8}\sin^2\theta-\frac9{1024}\sin^4\theta.
\]
Define \(\mathcal S_0=4\pi\int\rho^2(G\theta'^2/2+V_\alpha)d\rho=T+V\),
\(\mathcal I=4\pi\int\rho^2(3/4)\sin^2\theta\,d\rho,\ q=m_c\mathcal I\).
Stationarity implies
\[
 T+3V=0,\qquad \mathcal S_0=\frac23T>0,\qquad
 \mathcal S_0'(\alpha)=\frac{\mathcal I}{2}=\frac q{2m_c}.
\]
The leading binding \(B=\mathcal S_0-\alpha q/(2m_c)\) therefore obeys
\[
 B'(\alpha)=-\frac{\alpha}{2m_c}q'(\alpha).
\]
A charge fold is an extremum of binding, not necessarily its zero. The computed fold has positive B; binding crosses zero farther along the larger-charge branch.
The usual slope diagnostic \(dQ/d\omega<0\) translates to \(q'(\alpha)>0\) at leading order, but needs spectral assumptions. It is not substituted for curvature calculations.

For a smooth full-model branch, \(dE/dQ=\omega\), hence
\(d(E/Q)/dQ=-S_\omega/Q^2<0\). This alone does not exclude competing branches or prove nonlinear stability.
The reduced fold constants imply full-model minimum-charge scaling only conditionally on persistence of the branch/fold. A full-model fold-exponent sweep is still needed.

## Predicting the first correction

Follow the exact lambda-dependent valley:
\[
 f=\frac{\sqrt3}{2}\sin\theta,\qquad
 s=\sqrt{\frac3{16-\delta}}\cos\theta+\delta\eta.
\]
In the scaled action \(\mathcal S_\delta=\sqrt\delta(E-\omega Q)\), the order-delta density is
\[
 \frac3{512}\sin^2\theta\,\theta'^2-\frac9{16384}\sin^4\theta+
 s_0'\eta'+3\cos^2\theta\,\eta^2.
\]
This local expansion is checked symbolically. Radial integration by parts replaces \(s_0'\eta'\) with \(-(\Delta_\rho s_0)\eta\), assuming vanishing boundary terms. Quadratic minimization gives
\[
 \eta_*=\frac{\Delta_\rho s_0}{6\cos^2\theta},
\]
\[
 \mathcal S_1=4\pi\int\rho^2\left[
 \frac3{512}\sin^2\theta\,\theta'^2-\frac9{16384}\sin^4\theta
 -\frac{(\Delta_\rho s_0)^2}{12\cos^2\theta}\right]d\rho.
\]
Evaluation on the leading solution is sufficient: its first variation against the order-delta tangent shift vanishes. At fixed alpha,
\[
 \mathcal S_\delta=\mathcal S_0+\delta\mathcal S_1+O(\delta^2),
 \quad Q\delta^{3/2}=q_0+\delta q_1+O(\delta^2),
\]
\[
 q_1=2m_c\mathcal S_1'+\frac{3/64-\alpha}{m_c}\mathcal S_0'.
\]
The last term comes from \(Q\delta^{3/2}=2\omega\,d\mathcal S_\delta/d\alpha\) and must not be dropped. The coefficients come entirely from reduced solutions before new full-model evaluations, without fitting full-model outputs; parameter-difference refinement checks the derivatives.

This is a long-wavelength formal expansion, not a uniform theorem at either endpoint of the frequency window. The coordinate expansion is nonuniform as \(\cos\theta(0)\to0\). Its negative higher-gradient correction must not be promoted to a standalone ultraviolet theory.

## Correct fixed-charge curvature

The real-amplitude local block in physical radius is
\[
 L_+=-\nabla_l^2+
 \begin{pmatrix}
 1+3f^2+4s^2-\omega^2&8fs\\
 8fs&-3+4f^2+3\lambda s^2
 \end{pmatrix}.
\]
Only l=0 receives the additional
\[
 \frac{4\omega^2}{I}|(f,0)\rangle\langle(f,0)|
\]
with volume-weighted inner product. The charged-field phase block is
\(L_-=-\nabla_l^2+1+f^2+4s^2-\omega^2\).
There is no Sigma-phase block for real Sigma.
Our rho-coordinate matrices report physical curvature eigenvalues divided by delta. The transformation \(u=\rho\,\delta f\) symmetrizes the radial kinetic term; the rank-one vector uses the same discrete inner product.

Phase and translation zero modes are independent checks. Small negative symmetry-mode values must be refined before labeling them instabilities. A persistent negative radial mode establishes an energy saddle; positive sampled blocks do not prove nonlinear stability.

### Conditional analytic removal of nonradial energy instabilities

For an exact nodeless charged profile \(f>0\), \(L_- f=0\). The ground-state identity is
\[
 \int\eta L_-\eta=\int f^2|\nabla(\eta/f)|^2\ge0
\]
for admissible perturbations with vanishing boundary terms. Thus a nodeless solution has no negative phase-energy direction in this model, independently of small grid artifacts.

There is also a useful two-component identity. Suppose an exact spherical profile satisfies \(f,s>0,\ f'<0,\ s'>0\) for \(R>0\). Differentiating its equations gives \(L_{+,l=1}(f',s')=0\).
Flip the sign of the first component: the transformed operator has off-diagonal \(-C\), \(C=8fs>0\), and positive zero vector \(\psi=(-f',s')\). For perturbations \(\eta_i=\psi_i z_i\),
\[
 \langle\eta,\widetilde L_{+,1}\eta\rangle
 =4\pi\int R^2\left[
 \psi_1^2(z_1')^2+\psi_2^2(z_2')^2
 +C\psi_1\psi_2(z_1-z_2)^2\right]dR\ge0.
\]
This follows by radial integration by parts; the local divergence identity is checked symbolically in test_audit.py. The origin/infinity boundary terms must vanish, and strict interior positivity/monotonicity is an assumption, not a consequence of solver convergence.
For l>=2, add the positive diagonal centrifugal difference \([l(l+1)-2]/R^2\).
Accordingly, for this exact monotone class, any negative energy-curvature direction must lie in the fixed-charge l=0 amplitude sector. This is a conditional quadratic-energy result, not nonlinear or quantum stability and not a priority claim. It does not apply automatically to nodal, nonmonotone, gauged, or complex-Sigma extensions.

## Foundational parametric instability

For prescribed geometry, take \(\mu_1(A)=A,\mu_2=0,h_A=h_S=1,k_A=k_S=2\) and
\(\kappa^2(s)=1+\epsilon\cos\Omega s\).
The active internal Hessian is
\[
 K(s)=\begin{pmatrix}2&\kappa^2(s)\\\kappa^2(s)&2\end{pmatrix}.
\]
Its eigenvalues \(2\pm\kappa^2(s)\) remain positive for \(0<\epsilon<1\).
Yet \(x_\pm=(A\pm S)/\sqrt2\) obey
\[
 x_\pm''+[2\pm(1+\epsilon\cos\Omega s)]x_\pm=0.
\]
At \(\Omega=2\), the minus mode is resonant with the known Mathieu prediction \(\gamma=\epsilon/4+O(\epsilon^3)\).
The monodromy calculation includes detuned controls, determinant-one checks and a tighter integration check. It demonstrates the limitation of instantaneous-Hessian positivity, not a solution of the coupled worldline equations.

## Prior art and open bridges

Charge-stabilized solitons, constrained-energy analysis, effective-field elimination and parametric resonance are established methods.
Relevant primary sources include [Kim and Nugaev (2024)](https://link.springer.com/article/10.1140/epjc/s10052-024-13167-w), [Heeck and Sokhashvili (2023)](https://link.springer.com/article/10.1140/epjc/s10052-023-11710-9), and [Panin and Smolyakov (2017)](https://arxiv.org/abs/1612.00737). The latter specifically warns against exporting ordinary nongauged slope criteria to gauged theories.
Priority for this particular correction or reduced constants is unestablished; a focused quartic-extended-model comparison remains necessary.

The foundational constrained Hamiltonian, torsion boundary, quantization and observable map remain open. The two-field model has not been derived from that worldline action. These calculations establish no new element, planet, particle spectrum or uniquely SEFI mechanism.
