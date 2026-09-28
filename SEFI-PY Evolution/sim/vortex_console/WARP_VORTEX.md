# Warp-vortex addition

Opt-in NumPy module and CLI. Existing server, interface and computational functions are unchanged.

```powershell
python warp_vortex.py input.npz --output warp_result.npz
python -m unittest test_warp_vortex -v
```

Input arrays: `x,y,z` (increasing Cartesian axes, at least three nodes each); `D` (nx,ny,nz,3); `kappa,tau_g,envelope_level,sigma_weight` (nx,ny,nz); `sigma_residual` (nx,ny,nz,n_constraints). Optional `--identity '{"a":1}'` supplies scalar parameters only. Output paths must not already exist. NPZ loading disables pickle. This is a separate command alongside `server.py`, not a new browser tab.

- G_ij = partial_j D_i; H_ijk = partial_k partial_j D_i (component Hessian, not Laplacian).
- Omega_D = curl D; Omega_D,E = chi_E Omega_D.
- E = {e(D;Ia) <= 0}.
- Sigma membership: ||rSigma||_2 <= tolerance (default 1e-8).
- Effective weight = chi_Sigma w; w >= 0.
- Z_Sigma,E = integral chi_E chi_Sigma w |Omega_D|^2 dV.
- RMS_E = sqrt(integral chi_E |Omega_D|^2 dV / integral chi_E dV).

The caller supplies envelope levels, constraint residuals and weights; no framework-specific definitions are assumed. CLI arrays are pre-evaluated; Python callbacks can evaluate these quantities using D, its derivatives, invariants and Ia. Membership is not a dynamical stability result. Empty-envelope RMS is NaN in NPZ and null in the JSON summary. Quadrature uses tensor-product trapezoidal weights; derivatives use second-order finite differences with one-sided boundaries. No periodic or cylindrical differentiation is implied.

Worldline input uses arc length and a right-handed orthonormal Frenet frame with kappa > 0 and tau_g interpreted explicitly as Frenet torsion. Full Cartesian gradient samples are required: a line profile cannot determine transverse curl. Returned tangent and directional-derivative residuals expose input consistency. Existing axisymmetric console arrays must not be passed as Cartesian displacement data.

API: `displacement_components`, `displacement_derivatives`, `curl_displacement`, `compute_warp_envelope`, `check_stability_manifold`, `cartesian_volume_weights`, `stability_weighted_vortex`, `warp_vortex_status`, `worldline_warp_vortex_profile`.
