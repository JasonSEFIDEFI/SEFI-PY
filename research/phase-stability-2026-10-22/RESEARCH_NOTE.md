# From correction-versus-fatigue language to a measurable phase-stability model

Jason Duran Dutton, SEFI Research Program | October 22, 2026

**Status:** a proposed classical engineering extension with explicit added assumptions, not a derivation of a universal SEFI law. It does not replace submitted manuscripts. No biological, gravitational, consciousness or quantum-error-correction conclusion follows.

## 1. Source audit and corrections

`V:/V Drive and GitHub Content/Rough/1.tex` is itself the text-only transcription used to produce `1.pdf`; it is not an independent uncorrupted original. It explicitly preserves extracted ambiguities. No original source was edited. These corrections apply to its written definitions, not to an unverifiable earlier intended equation.

For Phi=A sin(omega t-kx+theta), Phi_x=-Ak cos(...) and Phi_xx=-Ak² sin(...). With the transcription's T=k Phi_x, R=T_x=-Ak³ sin(...), not the positive sign printed there. The exact curvature of the planar graph (x,Phi(x)) is |Phi_xx|/[1+Phi_x²]^(3/2); a signed convention omits the absolute value. The second derivative approximates signed graph curvature only in the small-slope regime and is not automatically spacetime curvature.

Small amplitude A tending to zero does not justify replacing sin(omega t-kx+theta) by its argument. That linearization requires a small argument modulo the chosen expansion point. A single sinusoid is a special solution, not the general solution of a wave equation or of every dynamic system.

If x is a physical length and k an inverse length, [T]=[Phi]/length², [R]=[Phi]/length³, while [integral |T| dx]=[Phi]/length. Therefore R>F as written needs dimensional coefficients or a nondimensional construction. A spatial integral of instantaneous load does not by itself establish accumulated irreversible damage or a failure-time prediction. A prescribed phase increment after threshold crossing implements a rule; it does not independently validate that rule.

## 2. Explicit new engineering assumptions

Let phi(t)=theta_signal(t)-theta_reference(t) be a measured, unwrapped **relative** phase; a fixed nominal carrier has already been removed. This is an optical/coherent-signal control variable, not an observable absolute phase of an unknown quantum state.

Let delta be constant free-running detuning in radians/second, K0 positive maximum correction gain, and d>=0 a dimensionless recoverable actuator-degradation state. Define

    K(d) = K0/(1+d),
    u = -K(d) sin(phi),
    dphi/dt = delta + u,
    dd/dt = alpha u² - beta d.

K0, delta and u have inverse-time units (radians treated as dimensionless); beta has inverse-time units; alpha has time units. Thus both differential equations are dimensionally consistent. Control activity accumulates degradation and beta describes recovery. At d=0, dd/dt>=0, so the physical half-plane is forward invariant for the continuous model. This is a phenomenological actuator-memory model, **not a deduction** that all material fatigue, aging or decoherence has this form.

Every new relation above is an added modeling choice. The reciprocal gain law and activity-squared loading need separate calibration and could be wrong. The constant-gain limit alpha=0 with d(0)=0 is the conventional Adler phase equation. Introducing d gives an ordinary nonlinear feedback system; it is not evidence of unique SEFI dynamics or improved control.

## 3. Derived stationary locking criterion

At a locked equilibrium, dphi/dt=0, so u*=-delta. Consequently

    d* = alpha delta² / beta,
    sin(phi*) = delta (1 + alpha delta²/beta) / K0.

The principal equilibrium exists and is locally asymptotically stable when

    |delta| (1 + alpha delta²/beta) < K0.

This is a precise replacement **within this model** for the vague comparison of correction capacity with accumulated stress. It is not the original dimensionally inconsistent R>F expression with relabeled variables.

Proof: put k*=K0/(1+d*) and c=cos(phi*)>0. The two-state Jacobian has

    trace = -k* c - beta - 2 alpha delta²/(1+d*) < 0,
    determinant = beta k* c > 0.

The characteristic polynomial therefore has roots with negative real parts. The alternate branch with negative cosine has negative determinant and is a saddle. At equality the determinant vanishes; the strict local-stability proof does not apply. Absence of a stationary equilibrium above the threshold does not by itself classify every possible nonlinear transient. The local result is not a global convergence theorem.

For alpha=0, this reduces to the familiar |delta|<K0. With K0=2/s, beta=1/s and alpha=1 s, stationary locking requires |delta|<1/s, whereas the ideal constant-gain model allows |delta|<2/s. The degradation extension predicts a *smaller* locking range, not superior control performance.

## 4. Nondimensional form

Set tau=beta t, Delta=delta/beta, k0=K0/beta and a=alpha beta. Then

    dphi/dtau = Delta - k0 sin(phi)/(1+d),
    dd/dtau = a [k0 sin(phi)/(1+d)]² - d.

The stationary condition becomes |Delta|(1+a Delta²)<k0. Every quantity in this comparison is dimensionless; the actuator and reference-frequency scales remain explicit when converting back to measurements.

## 5. A bounded experiment

Use a classical optical phase-locking apparatus or electronic phase-tracking testbed with measurable gain degradation and recovery. First measure K0 from the phase-detector/actuator response, beta from recovery after controlled loading, and alpha from loading transients. Fit those parameters on separate calibration records. Independently test whether gain really follows K0/(1+d); if not, reject or revise this ansatz before predicting the lock boundary.

Then sweep detuning on held-out runs and compare the predicted boundary with both constant-gain Adler dynamics and a conventional actuator-state model supplied the same measured gain/recovery information. Measure relative phase, applied control activity, gain and phase slips. Define a finite observation window, uncertainty and loss-of-lock criterion prospectively. Control for detector delay, saturation, temperature, gain variation and external noise. With stochastic noise, deterministic equilibrium existence alone does not guarantee indefinite locking; noise-driven slips need a separate analysis.

A physical d estimate can be tied to independent gain measurements by d=K0/K_measured-1 within the gain law's valid regime. That identification does not license fitting d from the very phase failures being predicted. An internal software d state without an independent optical/actuator meaning would only test the simulator.

If this law predicts failures better than an incorrectly idealized constant-gain model, that supports including actuator memory. It does not establish superiority to an equivalent conventional degradation model or confirm the full SEFI ontology. A laser-only test establishes classical phase tracking, not retention of unknown photonic qubits.

## 6. Existing work and novelty boundary

Noisy phase locking and phase slips already have established models. Relevant primary sources include Freund et al., *Noise-Induced Phase Synchronization: Theoretical and Experimental Results*, https://toniq.c2n.universite-paris-saclay.fr/media/filer_public/9e/a9/9ea92df4-d547-427a-8814-cca9adb5acec/freundfnl03.pdf ; and the oscillator study *Phase-locked spiking of inner ear hair cells and the driven noisy Adler equation*, https://pmc.ncbi.nlm.nih.gov/articles/PMC4213446/ . Optical phase-lock cycle slips are discussed in the experimental optical-clock work hosted by NIST, https://tf.nist.gov/general/pdf/2206.pdf . These references establish overlap, not priority or novelty for the particular actuator ansatz. A systematic novelty review has not been completed.

The contribution of this note is making one proposed mechanism explicit and testable, documenting its assumptions and analytical consequences. No new physical law, superiority, hardware validation, or journal acceptance is claimed.
