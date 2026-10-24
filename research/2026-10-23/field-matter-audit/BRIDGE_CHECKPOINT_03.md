# Five-minute checkpoint 03: internal charge to spacetime current

23 October 2026. This constructs a distributional spacetime current from a
worldline/internal charge. It is a conditional mathematical bridge, not a
smooth bulk field, a quantum amplitude or a matter-production model. No
existing action is changed.

## 1. Definition and assumptions

Let y(s) be a smooth future-directed timelike worldline, with proper time s,
tangent T^mu=dy^mu/ds, and a real continuously differentiable internal charge
J(s). J may be either candidate charge of checkpoint 02, with its actual
proper-time derivative supplied by the corresponding torque identity. For a
finite interval [a,b], define

    j^mu(x) = integral_a^b ds J(s) T^mu(s) delta^4(x-y(s)).

This means, for every smooth compactly supported covector test field A_mu,

    <j,A> = integral_a^b ds J(s) T^mu(s) A_mu(y(s)).

The delta is the spacetime coordinate delta, not a probability law. With
dimensionless J, this current has dimensions length^-3 in natural units.
For an orientation-preserving parameter u the same expression is

    <j,A> = integral du J(s(u)) (dy^mu/du) A_mu(y(u)).

Thus proper-time factors cancel correctly; J is composed with the parameter
change and must not be assigned an extra density factor. This is a covariant
worldline vector distribution in the assumed flat spacetime. It does not
derive spacetime geometry. Numerical tests use coordinate time as parameter
on timelike curves; this is justified by the line-integral expression above.

## 2. Derive the divergence, including its signs

For a smooth compactly supported scalar test f, the definition of a
distributional derivative and the chain rule give

    <partial_mu j^mu,f>
      = -integral_a^b ds J T^mu partial_mu f(y)
      = -integral_a^b ds J d[f(y)]/ds
      = J(a) f(y(a)) - J(b) f(y(b))
        + integral_a^b ds J_dot(s) f(y(s)).

Consequently

    partial_mu j^mu(x)
      = J(a) delta^4(x-y(a)) - J(b) delta^4(x-y(b))
        + integral_a^b ds J_dot(s) delta^4(x-y(s)).

The start sign is positive and the end sign negative. Substituting checkpoint
02's torque formulas gives the source explicitly; no new conservation law is
assumed. On an isolated smooth curve away from endpoints, J_dot=0 implies
local conservation. Conversely, testing against arbitrary functions near the
curve detects any nonzero J_dot. Other sectors or overlapping currents could
cancel sources, but such sectors must be specified rather than presumed.

Even a constant nonzero charge on a **finite isolated segment** has endpoint
sources. For complete histories with no endpoints in the support of f, or
where endpoints balance with other explicitly included currents, these terms
vanish. Truncating an otherwise conserved history to an observation window
creates bookkeeping boundary terms; it is not evidence of physical charge
creation or failure of the underlying conservation law.

## 3. Charge on an inertial time slice

Because the worldline is future-directed timelike, y^0(s) is strictly
increasing and crosses a fixed inertial time slice at most once. If it crosses
at s_t in the interior of the segment, integrating j^0 over that slice gives

    Q(t) = integral d^3x j^0(t,x)
         = integral ds J(s) (dy^0/ds) delta(t-y^0(s))
         = J(s_t).

The positive delta-function Jacobian cancels dy^0/ds. There is no extra
Lorentz-factor multiplier for the slice charge. This statement assumes the
entire spatial slice is included; a bounded spatial region also has transport
flux through its boundary. The slice-charge identity is derived here; no
numerical delta-kernel approximation is claimed as a separate test.

## 4. Endpoint networks: what can balance and what is not established

For a collection of segments, the coefficient at a shared event is

    sum_outgoing J_start - sum_incoming J_end.

It vanishes precisely when incoming and outgoing charge totals match. Exact
rational checks verify that incoming charge 3 splitting into outgoing charges
1 and 2 balances, while splitting into 1 and 1 leaves coefficient -1.
Two outgoing segments with charges +1 and -1 also cancel their common initial
endpoint coefficient. Their remote endpoints still need accounting.

This is an algebraic consistency condition on a specified collection of
currents. A branching network is additional structure beyond a single smooth
worldline; no branching dynamics is derived. Cancellation for opposite signs
does not demonstrate pair production, identify antiparticles, provide energy
or establish a transition amplitude. No future-directed timelike curve is
silently reversed in time.

## 5. Actual independent checks

Four new tests cover:

1. Endpoint signs against the analytic polynomial result -8.16 on the finite
   segment y(t)=(t,0.2t,0.1t,0), 0<=t<=2, J=2, f=t^2+x*y.
2. A test vanishing at both endpoints, f=t^2(2-t)^2: constant charge gives zero
   weak interior divergence; J(t)=1+t/2 gives exactly 8/15 from its interior
   source. Numerical results were -7.94e-17 and 0.5333333333333333.
3. A nonlinear orientation-preserving reparameterization t=u+0.2 sin(u) on a
   curved timelike trajectory: weak current pairings agree to 1e-11, including
   a scalar-gradient test of the weak divergence.
4. Exact endpoint charge balance for the segment networks described above.

The polynomial tests can be multiplied by a smooth compact cutoff equal to
one on a neighborhood of the finite curve, leaving the evaluated values and
gradients unchanged. Thus the computations test legitimate distributional
identities without claiming that a bare polynomial is compactly supported.

**All 17 local tests passed in 4.25 seconds**, including the preceding 13.
Quadrature error estimates and actual residuals are saved in
`worldline-current-results.json`; the constant-charge residual is at roundoff
scale, not measured charge violation. These tests concern weak distributions,
not a gridded charge-deposition scheme or a PDE solver. The existing CI
workflow collects the new tests; remote success is not inferred from this run.

Reproduce from repository root:

```sh
python -m pytest research/2026-10-23/field-matter-audit -q
python research/2026-10-23/field-matter-audit/worldline_current.py
```

The science export refuses to overwrite its existing result file.

## 6. Consequence for the research program

We now have a defined connection between an internally conserved quantity and
a locally conserved spacetime current, under explicit boundary assumptions.
The current remains singular and supported on the supplied trajectory. It
does not yield the scalar soliton equations, electromagnetic charge,
fermionic spin, quantum statistics, a field-generation mechanism or a smooth
matter distribution. Conservation restricts possible couplings; it does not
choose a coupling or predict a particle.

Next five-minute task: examine a proposed vector-field source coupling
integral j^mu A_mu and derive its response to A_mu -> A_mu+partial_mu chi.
Use the endpoint and torque terms above to determine when that proposed
coupling is gauge-invariant and what additional source accounting it would
require. This is a consistency test of a clearly labeled possible extension,
not an identification of J with electric charge.
