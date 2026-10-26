# One-session continuation: full-model charge folds

26 October 2026. Dimensionless classical two-field model; not a particle or nuclear model derived from the foundational SEFI action.

## Prediction before full-model evaluation

The previous audit derived
\[
 q(\alpha,\delta)=Q\delta^{3/2}
 =q_0(\alpha)+\delta q_1(\alpha)+O(\delta^2),
\quad \delta=16-\lambda,
\]
with \(q_0=2m_c\mathcal S_0'\) and
\[
 q_1=2m_c\mathcal S_1'+
 \frac{3/64-\alpha}{m_c}\mathcal S_0',\quad m_c=\sqrt{7/4}.
\]

Let the reduced fold be \(\alpha_0\), with \(q_0'(\alpha_0)=0\) and \(q_0''(\alpha_0)>0\).
A Taylor expansion of the stationarity condition gives
\[
 \alpha_{\mathrm f}=\alpha_0+\delta\alpha_1+O(\delta^2),\qquad
 \alpha_1=-q_1'(\alpha_0)/q_0''(\alpha_0),
\]
\[
 Q_{\min}(\delta)=\delta^{-3/2}
 [q_0(\alpha_0)+\delta q_1(\alpha_0)+O(\delta^2)].
\]
The shift of alpha does not contribute to the charge minimum at first order because \(q_0'(\alpha_0)=0\). At the fold,
\[
 q_1'=2m_c\mathcal S_1''-\mathcal S_0'/m_c,
\]
since \(\mathcal S_0''=q_0'/(2m_c)=0\).
Three parameter-difference spacings check derivative sensitivity. Coefficients are saved before the new full-model evaluations, without fitting to those evaluations. This is not an external preregistration.

## Full-model boundary-value problem

The original potential is
\[
 U=f^4/4+f^2/2+2f^2s^2+\lambda(s^2-3/\lambda)^2/4.
\]
With \(\rho=\sqrt\delta R,\ \omega^2=m_\lambda^2-\delta\alpha\),
\[
 f''+2f'/\rho=(1-\omega^2+f^2+4s^2)f/\delta,
\]
\[
 s''+2s'/\rho=(-3+4f^2+\lambda s^2)s/\delta.
\]
Boundary data: \(f'(0)=s'(0)=0,\ f(L)=0,\ s(L)=\sqrt{3/\lambda}\).
The additional central-amplitude condition \(f(0)=(\sqrt3/2)\sin c\) determines the unknown alpha.
Minimize the computed scaled charge over c on the selected positive monotone spherical branch. This is not minimization over every possible stationary or nonspherical solution.

Check residuals, an independent off-node ODE residual, frequency admissibility, monotonicity, the virial identity, and positive finite-difference curvature of charge at the minimum. A longer-domain/tighter-tolerance repeat checks one fold.

## Why a charge fold has a radial energy zero mode

For a smooth stationary family \(F(c)=(f(c),s(c))\), the fixed-frequency Hessian obeys
\[
 H_\omega F_c=2\omega\omega_c(f,0).
\]
The fixed-charge Hessian adds \(4\omega^2|(f,0)\rangle\langle(f,0)|/I\), where \(I=\int f^2\). Consequently
\[
 H_QF_c=
 2\omega\left(\omega_c+\omega I_c/I\right)(f,0)
 =\frac{2\omega Q_c}{I}(f,0).
\]
Thus, at a charge fold \(Q_c=0\), a nonzero family tangent is a null direction of the radial fixed-charge Hessian. The assertion requires smoothness, admissible decay and a nonzero tangent. A simple index crossing additionally requires nondegeneracy; the identity alone does not establish it.

In scaled coordinates the curvature operator is \(H_Q/\delta\), and
\[
 (H_Q/\delta)F_c=\frac{2\omega q_c}{\delta I_\rho}(f,0),
\quad I_\rho=4\pi\int\rho^2 f^2d\rho.
\]
The factors of \(4\pi\) cancel in \(q_c/I_\rho\). The radial matrix's rank-one normalization uses its matching discrete inner product.

The session tests radial eigenvalues below, at and above the charge fold, at two couplings and two spatial resolutions. At the fold it also compares the lowest eigenvector to the numerically differentiated stationary-family tangent. Small finite-grid negative eigenvalues at a zero mode must be refined before interpretation.

### Independent family derivative

For \(u=f_c,v=s_c\) and the unknown parameter \(\alpha_c\), differentiating the scaled equations gives
\[
 u''+2u'/\rho=A u+C v+\alpha_c f,\qquad
 v''+2v'/\rho=C u+D v,
\]
where \(A=(1-\omega^2+3f^2+4s^2)/\delta,\ C=8fs/\delta,\ D=(-3+4f^2+3\lambda s^2)/\delta\).
The differentiated boundary data are \(u(0)=(\sqrt3/2)\cos c,\ u'(0)=v'(0)=0,\ u(L)=v(L)=0\).
Then \(\omega_c=-\delta\alpha_c/(2\omega)\),
\(I_{\rho,c}=8\pi\int\rho^2fu\,d\rho\), and \(q_c=\omega_c I_\rho+\omega I_{\rho,c}\).

This linear BVP avoids subtracting nearby nonlinear profiles. It also supports a derivative-based correction to the minimizing central parameter.
The family identity must be checked with matching outer boundary conditions. Shortening the operator domain while keeping an unmodified longer-domain tangent introduces a boundary residual that can grow under grid refinement. The original mismatched diagnostics are retained and the matching-domain certificate uses three resolutions.

## Limits

This follows ordinary constrained-energy/soliton methods. Novelty is unestablished. Positive quadratic energy curvature is not a nonlinear stability theorem, and a charge minimum does not mark the free-particle energetic-binding threshold.
The model's free-field mass, charge and coupling are dimensionless model quantities. No measured matter, new element, physical spectrum, or unique SEFI advantage is inferred.
