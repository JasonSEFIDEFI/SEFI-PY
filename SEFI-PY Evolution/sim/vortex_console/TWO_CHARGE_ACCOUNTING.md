# Two-charge initial-data test — 24 September 2026

The fixed-carrier-charge torus has a large original-field charge deficit. We constructed initial velocities that restore that charge while retaining carrier charge Q=1000. This is a charge-matched initial-data test, not a stationary two-charge solution, a physical evolution, or a stability proof.

Let F be the original-field profile, f a compact compensation window, and Delta Q_Phi its initial excess charge. Add

    d_T F = i deltaOmega f F,
    deltaOmega = -Delta Q_Phi / A,
    A = integral f |F|^2,
    D = integral f^2 |F|^2,
    chi = A^2/D.

The compensated physical background-subtracted energy is exactly

    E_matched = R_Q + (Delta Q_Phi)^2/(2 chi).

We used f=max(1-(r^2+z^2)/R^2,0)^4. The refined numerical spatial fields were extended by their homogeneous outer values. The grid spacing was .4. Both charge identities and the direct versus reconstructed energy were checked to absolute tolerance 1e-8.

| Compensation support R | Extra energy cost | E_matched/Q |
|---:|---:|---:|
| 20 | 3139.2354 | 4.047159 |
| 28 | 876.0322 | 1.783956 |
| 36 | 394.1885 | 1.302112 |
| 48 | 165.1574 | 1.073081 |
| 64 | 69.9013 | .977825 |

The carrier charge is 1000 in all cases, and the remaining original-field excess charge is at most 2.1e-12 in magnitude. A below-threshold comparison is possible in the largest tested surrounding compensation region. This does not establish an intrinsic localized mass: the prepared state includes that surrounding phase-velocity distribution. Axial momentum matching, fully constrained stationary solutions, non-axisymmetric disturbances and physical evolution remain unresolved. Axial angular momentum is unaffected by the axisymmetric correction to the original field's phase velocity.

This identifies a precise next question: can a localized excitation and its background be treated self-consistently at matching conserved charges, with stable, reproducible energy and momentum response? Do not present the current energy comparison as an observed particle or as proof of dynamical binding.
