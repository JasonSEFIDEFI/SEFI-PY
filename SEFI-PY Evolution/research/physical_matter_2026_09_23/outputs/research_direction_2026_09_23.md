# Toward physical matter from field excitations

Jason Duran Dutton | Research direction and mathematical checkpoint | 23 September 2026

## Explicit goal and present conclusion

The explicit goal is to derive physical matter from a mathematically consistent field framework motivated by a primordial unity hypothesis, and to determine whether its excitations can generate the spacetime geometry universally experienced by matter. The primordial unity and its proposed Big Bang transformation are motivating hypotheses. They are not established premises of the calculations below.

We have derived an effective metric for one low-energy phase sector and constructed numerical moving vortex rings, including a provisional two-field charged loop. We have not derived an observed particle, universal gravitational dynamics, quantum statistics, or a cosmological history. The latest work provides limited dynamical evidence for the moving loop and an analytical obstruction to a particular rest-state interpretation. This is a research checkpoint, not a submission-ready claim of physical matter.

The earlier SEFI, DEFI and GWFM manuscripts supply the research motivation and questions about identity, coupling and worldlines. They are not used here as independently established physical laws. Their proposed structures require separate variational and empirical checks before inclusion in this model.

## 1 Starting field theory

We assume flat background spacetime with signature (+---), set its limiting speed to one, and use a canonical complex scalar field. With positive parameters alpha and beta,

$$
\mathcal L_\Phi=\frac12\partial_\mu\Phi^*\partial^\mu\Phi-\frac\alpha2|\Phi|^2-\frac\beta4|\Phi|^4.
$$

The field equation and homogeneous rotating solution are

$$
\Box\Phi+(\alpha+\beta|\Phi|^2)\Phi=0,\qquad \Phi_\infty=R_\infty e^{i\Omega t},\qquad R_\infty^2=\frac{\Omega^2-\alpha}{\beta}.
$$

This construction assumes a spacetime and a rotating background. It does not derive either from a pre-spacetime entity. All numerical quantities below are dimensionless; no measured mass or length scale has been calibrated.

Writing the field in amplitude and phase variables and eliminating slow amplitude response gives the infrared phase theory, on the branch X greater than alpha,

$$
X=\partial_\mu\theta\partial^\mu\theta,\qquad P(X)=\frac{(X-\alpha)^2}{4\beta},\qquad K^{\mu\nu}=\frac{(X-\alpha)\eta^{\mu\nu}+2w^\mu w^\nu}{\beta},\quad w_\mu=\partial_\mu\theta.
$$

The linear phase equation is a wave equation with this effective kinetic tensor. One compatible metric in four dimensions is

$$
g_{\mu\nu}=\frac{\sqrt{(X-\alpha)(3X-\alpha)}}{\beta}\left(\eta_{\mu\nu}-\frac{2w_\mu w_\nu}{3X-\alpha}\right).
$$

For alpha=beta=1 and Omega=sqrt(2), the phase sound speed squared is 1/5. This is a sector-dependent effective geometry. The amplitude mode has a different low-momentum dispersion, and no Einstein equation or common metric for all excitations has been obtained. The phase-only approximation must not be extended through a vortex core where its elimination assumptions fail. This distinction between effective propagation geometry and gravitational dynamics is central to analogue-gravity research [1].

## 2 Moving vortex rings

For the single-field traveling ansatz, a rescaled longitudinal coordinate reduces the profile equation exactly to a Gross-Pitaevskii traveling-wave equation. In the numerical normalization used here,

$$
Z=\gamma(z-vt),\quad \gamma=(1-v^2)^{-1/2},\quad \Phi=e^{i\Omega t}\psi(r,Z),\quad \nabla^2\psi+ic\partial_Z\psi+(1-|\psi|^2)\psi=0,\quad c=2\Omega\gamma v.
$$

Known Gross-Pitaevskii vortex-ring existence results motivate this profile construction [2]. Their time-dependent stability results cannot be imported into this second-order relativistic field equation. Earlier finite-domain numerical rings showed winding one and restricted persistence tests; restricted fixed-momentum and fixed-charge trial variations did not constitute full stability proofs.

## 3 Explicit two-field extension

To investigate localized charge, a second canonical complex field Sigma was added as a hypothesis. It is not a derivation of primordial duality. The potential is

$$
V=\frac12|\Phi|^2+\frac14|\Phi|^4+\frac\lambda4|\Sigma|^4+\frac12(4|\Phi|^2-3)|\Sigma|^2.
$$

The baseline uses lambda=25. The carrier has bulk mass squared 1 and negative local mass squared -3 where the original amplitude vanishes. A negative eigenvalue of the uncondensed straight-core carrier operator motivated a nonlinear coupled calculation. At carrier frequency 0.3, a refined straight-core solution has carrier amplitude 0.310063 at its center and charge per unit length 0.282372. Positive radial real-amplitude Hessian eigenvalues cover only that restricted variation, not full dynamical stability.

Charge-carrying string cores have theoretical precedent [3]. Vorton existence and stability in other models do not establish these properties in this ungauged, rotating-background model [4].

For the closed moving loop, use

$$
\Phi=e^{i\Omega t}\psi(r,Z),\qquad \Sigma=e^{i\nu\gamma(t-vz)}s(r,Z),\qquad Z=\gamma(z-vt).
$$

With real s, the exact coupled profile equations are

$$
\nabla^2\psi+ic\partial_Z\psi+(1-|\psi|^2-4s^2)\psi=0,\qquad \nabla^2s+(3+\nu^2-4|\psi|^2-\lambda s^2)s=0.
$$

The conserved carrier charge is Q = nu times the integral of s squared over the comoving spatial coordinates. Exponential carrier localization requires absolute nu less than 1. The charge is a continuous classical global Noether charge; it has not been identified with electric charge.

At c=0.6, nu=0.3 and lambda=25, speed v=0.207514. The following simultaneous two-field solutions were obtained on finite cylinders with radial extent and axial half-length L:

| L | Grid spacing | Ring radius | Carrier charge | Maximum discrete residual |
|---|---|---|---|---|
| 16 | 0.666667 | 5.366207 | 7.085225 | 1.3e-12 |
| 16 | 0.500000 | 5.357405 | 7.133380 | 3.9e-11 |
| 24 | 0.666667 | 5.157070 | 6.920481 | 6.6e-13 |

Each profile has meridional winding +1. The larger-domain radius differs by approximately 3.9%, so infinite-domain convergence is not established. A small discrete residual verifies equation solving on a grid, not physical accuracy or stability.

## 4 New time-dependent checks

The full complex two-field equations were evolved in moving coordinates, including both phase and amplitude response. No field component was held fixed. The unperturbed control remains within approximately 6e-9 in the reported weighted field norms through time 100.

A Gaussian real perturbation of peak 0.001 was added to both fields, with initially zero moving-frame velocities. This slightly changes the initial total carrier charge, which is then conserved within each run. At time 100 the loop retains winding one. Halving the time step from 0.025 to 0.0125 changes the final fields by at most 1.9e-8 pointwise; the refined run's maximum relative carrier-charge drift is 5.2e-13.

The raw carrier-field distance grows from 0.01028 to 0.43265. Most of this is global phase drift of approximately 0.08872 radians: after phase alignment, the final distance is 0.00870, and the carrier-amplitude distance is 0.00739. These results are consistent with small deformation in this particular test. They do not prove stability or a lifetime.

Two linear initial-value probes of non-axisymmetric sectors were also run. For azimuthal mode 1, the field-norm ratio peaks at 1.7725 and ends at 1.1083; for mode 2 it peaks at 1.2085 and ends at 1.0794. These sample one initial condition per sector, use reflecting finite boundaries and a fixed spatial grid, and are neither full eigenvalue searches nor nonlinear three-dimensional simulations. They do not exclude unexcited unstable modes or slow growth.

## 5 A rest-state obstruction

The strongest new analytical result concerns profiles at rest relative to the homogeneous background. Let Phi=exp(i Omega t) F(x) and Sigma=exp(i nu t) S(x), with canonical gradients, localized deviations and finite background-subtracted energy. Set x=|F| squared and y=|S| squared. For lambda=25, the frequency-constrained potential is exactly

$$
U_\nu(x,y)=\frac14(x-1+4y)^2+\frac94y^2+\frac{1-\nu^2}{2}y.
$$

For absolute nu at most 1 this is nonnegative. The profile functional is F=K+integral U, with K the nonnegative spatial gradient energy. Under a uniform dilation of a localized three-dimensional configuration,

$$
\mathcal F(\ell)=\ell K+\ell^3\mathcal U,\qquad \left.\frac{d\mathcal F}{d\ell}\right|_{\ell=1}=K+3\mathcal U=0.
$$

Nonnegativity forces the homogeneous configuration. This is a Derrick-type scaling obstruction [5], applicable to critical points of the frequency-constrained functional, including fixed-charge extrema with the frequencies as Lagrange multipliers. Spatial phase winding contributes to positive gradient energy and does not remove it.

The result excludes a nontrivial localized time-harmonic rest profile for these parameters and assumptions. It does not exclude the computed moving loop, an infinite straight string, more general time dependence, different interactions, or a rest frame with a flowing asymptotic background. Boosting a moving loop changes that background, so there is no contradiction. It is not a proof that all possible matter models based on fields fail.

## 6 A justified parameter test and its outcome

Keeping the same fields and changing the carrier self-coupling is an explicit new model choice. There is an analytically accessible window

$$
9<\lambda<16,\qquad \sqrt{\lambda}-3<\nu^2<1.
$$

In this window the zero-carrier-frequency bulk state remains a global minimum of its rotating-frame potential, while the frequency-constrained potential can become negative. To see the latter, set x=0 and y=(3+nu squared)/lambda: the potential equals one quarter minus (3+nu squared) squared divided by 4 lambda. Thus this particular positivity obstruction is removed without making the far-field carrier mass squared negative.

At lambda=12 and nu=0.85, the trial potential is -0.0386876. Two initial searches for rest loops with carrier azimuthal winding were attempted. The winding-1 search converged to the homogeneous state with negligible carrier, not a loop. The winding-3 search stalled with residual about 0.689. Neither establishes a rest-loop solution, and neither proves nonexistence in the revised model. This is the stopping point; further parameter or continuation searches should follow a method review.

## 7 What would count as deriving physical matter

The long-term goal requires more than a localized classical field configuration. A defensible route must establish the relevant solution branch and stability, derive its energy and momentum response, specify quantization and spin/statistics, identify measured charges and interactions, and produce testable scales or dimensionless predictions. A proposed universal geometry must govern all relevant excitations consistently and supply tested gravitational dynamics. A primordial origin claim additionally needs a cosmological evolution and comparison with observations.

The immediate research decision is whether to develop the revised two-field rest branch with constrained continuation or to retain the moving-ring model as an effective-medium study. Either route should first address full stability and the mismatch between sector-dependent propagation geometries. No result currently licenses the statement that these excitations produce the spacetime universally experienced by physical matter.

## 8 Computational implementation and provenance

The numerical work uses independent NumPy research scripts: nonlinear finite-difference Newton solvers, explicit fourth-order time stepping, and limited linear perturbation probes. A separate opt-in research directory preserves profiles, numerical histories, equations and failed searches. It does not alter the SEFI-PY engine's modules, dependencies, configuration or entry points. It is research software, not a validated extension of the engine's physical claims.

Calculations and this synthesis were developed with AI assistance. The included derivations and numerical evidence require independent expert review. The reproducibility archive distinguishes analytical conclusions, finite-grid results and unsuccessful searches; no external replication or journal acceptance is claimed.

## References

[1] C. Barcelo, S. Liberati and M. Visser, Analogue Gravity, Living Reviews in Relativity 14, 3 (2011). https://doi.org/10.12942/lrr-2011-3

[2] F. Bethuel, G. Orlandi and D. Smets, Vortex rings for the Gross-Pitaevskii equation, Journal of the European Mathematical Society 6, 17-94 (2004). https://ems.press/journals/jems/articles/92

[3] E. Witten, Superconducting strings, Nuclear Physics B 249, 557-592 (1985). https://doi.org/10.1016/0550-3213(85)90022-7

[4] R. A. Battye, S. J. Cotterill and J. A. Pearson, A detailed study of the stability of vortons, Journal of High Energy Physics 04, 005 (2022). https://arxiv.org/abs/2112.08066

[5] G. H. Derrick, Comments on Nonlinear Wave Equations as Models for Elementary Particles, Journal of Mathematical Physics 5, 1252-1254 (1964). https://doi.org/10.1063/1.1704233
