"""Opt-in Cartesian warp/curl diagnostics. Existing console modules are untouched."""
import argparse
import json
from pathlib import Path
import numpy as np


def _finite(value, name):
    a = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(a)):
        raise ValueError(name + ' must be finite')
    return a


def _parameters(Ia):
    for key, value in Ia.items():
        if not isinstance(key, str) or np.ndim(value) != 0:
            raise ValueError('Ia must map names to finite scalar parameters')
        _finite(value, 'Ia')


def _axis(value, name):
    a = _finite(value, name)
    if a.ndim != 1 or len(a) < 3 or np.any(np.diff(a) <= 0):
        raise ValueError(name + ' needs at least 3 increasing coordinates')
    return a


def _field(D, coordinates):
    if len(coordinates) != 3:
        raise ValueError('Three Cartesian axes required')
    axes = tuple(_axis(a, name) for a, name in zip(coordinates, ('x', 'y', 'z')))
    D = _finite(D, 'D')
    if D.shape != tuple(len(a) for a in axes) + (3,):
        raise ValueError('D shape must be (nx, ny, nz, 3)')
    return D, axes


def displacement_components(D, T, N, B):
    """Dt=D.T; Dn=D.N; Db=D.B; basis: right-handed orthonormal."""
    D = _finite(D, 'D')
    basis = np.stack([_finite(a, 'basis') for a in (T, N, B)], axis=-1)
    if D.shape[-1:] != (3,) or basis.shape != D.shape + (3,):
        raise ValueError('Basis and D shapes must agree')
    if not np.allclose(np.swapaxes(basis, -1, -2) @ basis, np.eye(3), atol=1e-7, rtol=0):
        raise ValueError('Basis must be orthonormal')
    if not np.allclose(np.linalg.det(basis), 1, atol=1e-7, rtol=0):
        raise ValueError('Basis must be right-handed')
    return tuple(np.sum(D * basis[..., i], axis=-1) for i in range(3))


def displacement_derivatives(D, coordinates):
    """G[...,i,j]=partial_j D_i; H[...,i,j,k]=partial_k partial_j D_i."""
    D, axes = _field(D, coordinates)
    G = np.stack(np.gradient(D, *axes, axis=(0, 1, 2), edge_order=2), axis=-1)
    H = np.stack(np.gradient(G, *axes, axis=(0, 1, 2), edge_order=2), axis=-1)
    return G, H


def curl_displacement(G):
    """Omega_D_i=epsilon_ijk G_kj."""
    G = _finite(G, 'gradient')
    if G.shape[-2:] != (3, 3):
        raise ValueError('Gradient must end in (3,3)')
    return np.stack((G[..., 2, 1]-G[..., 1, 2],
                     G[..., 0, 2]-G[..., 2, 0],
                     G[..., 1, 0]-G[..., 0, 1]), axis=-1)


def compute_warp_envelope(D, Ia, envelope_level):
    """E(D)={x:e(D;Ia)<=0}; e must return one value per grid point."""
    D = _finite(D, 'D')
    if D.shape[-1:] != (3,):
        raise ValueError('D must end in 3')
    _parameters(Ia)
    level = _finite(envelope_level(D, Ia), 'envelope level')
    if level.shape != D.shape[:-1]:
        raise ValueError('Envelope shape mismatch')
    return level <= 0


def check_stability_manifold(D, G, H, kappa, tau_g, Ia, sigma_residual, tolerance):
    """Sigma: ||rSigma(D,G,H,kappa,tau_g;Ia)||_2<=tolerance."""
    # Membership test, not a dynamical stability test.
    D = _finite(D, 'D'); G = _finite(G, 'G'); H = _finite(H, 'H')
    kappa = _finite(kappa, 'kappa'); tau_g = _finite(tau_g, 'tau_g')
    _parameters(Ia)
    shape = D.shape[:-1]
    if D.shape[-1:] != (3,) or G.shape != D.shape+(3,) or H.shape != D.shape+(3,3):
        raise ValueError('Derivative shape mismatch')
    if kappa.shape != shape or tau_g.shape != shape or np.any(kappa < 0):
        raise ValueError('Invariant grid shapes must agree; kappa >= 0')
    if not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError('Tolerance must be finite and nonnegative')
    r = _finite(sigma_residual(D, G, H, kappa, tau_g, Ia), 'Sigma residual')
    if r.ndim != len(shape)+1 or r.shape[:-1] != shape or r.shape[-1] < 1:
        raise ValueError('Sigma residual must have shape grid+(n_constraints,)')
    return np.linalg.norm(r, axis=-1) <= tolerance


def cartesian_volume_weights(coordinates):
    """Tensor-product trapezoidal integration weights."""
    if len(coordinates) != 3:
        raise ValueError('Three Cartesian axes required')
    weights = []
    for a in coordinates:
        a = _axis(a, 'axis'); d = np.diff(a)
        weights.append(np.r_[d[0]/2, (d[:-1]+d[1:])/2, d[-1]/2])
    return weights[0][:, None, None]*weights[1][None, :, None]*weights[2][None, None, :]


def stability_weighted_vortex(omega_D, sigma_weight, envelope_mask, dV):
    """Z_Sigma,E=sum(chi_E w_Sigma |Omega_D|^2 dV)."""
    omega_D = _finite(omega_D, 'Omega_D'); shape = omega_D.shape[:-1]
    w = _finite(sigma_weight, 'Sigma weight'); dV = _finite(dV, 'dV')
    mask = np.asarray(envelope_mask)
    if omega_D.shape[-1:] != (3,) or w.shape != shape or dV.shape != shape or mask.shape != shape:
        raise ValueError('Diagnostic shape mismatch')
    if mask.dtype != bool or np.any(w < 0) or np.any(dV <= 0):
        raise ValueError('Boolean mask, nonnegative weights and positive dV required')
    return float(np.sum(mask*w*np.sum(omega_D**2, axis=-1)*dV))


def warp_vortex_status(D, coordinates, kappa, tau_g, Ia,
                       sigma_residual, sigma_weight, envelope_level, tolerance=1e-8):
    """Omega_D=curl D; Omega_D,E=chi_E Omega_D; Z_Sigma,E=integral_E w_Sigma |Omega_D|^2."""
    D, axes = _field(D, coordinates)
    G, H = displacement_derivatives(D, axes)
    omega = curl_displacement(G)
    inside = compute_warp_envelope(D, Ia, envelope_level)
    sigma = check_stability_manifold(D, G, H, kappa, tau_g, Ia, sigma_residual, tolerance)
    w = _finite(sigma_weight(D, G, H, kappa, tau_g, Ia), 'Sigma weight')
    if w.shape != inside.shape or np.any(w < 0):
        raise ValueError('Sigma weights must match grid and be nonnegative')
    w = w*sigma
    dV = cartesian_volume_weights(axes)
    volume = float(np.sum(inside*dV))
    integral = stability_weighted_vortex(omega, np.ones_like(w), inside, dV)
    return dict(gradient=G, hessian=H, omega_D=omega, omega_envelope=inside[..., None]*omega,
                envelope_mask=inside, sigma_mask=sigma, effective_sigma_weight=w,
                envelope_volume=volume, rms_envelope=float(np.sqrt(integral/volume)) if volume else float('nan'),
                weighted_vortex=stability_weighted_vortex(omega, w, inside, dV))


def worldline_warp_vortex_profile(s, x, Dt, Dn, Db, T, N, B,
                                  kappa, tau_g, Ia, gradient_at_worldline):
    """D'= (Dt'-kappa Dn)T+(Dn'+kappa Dt-tau_g Db)N+(Db'+tau_g Dn)B."""
    # s: arc length; tau_g: Frenet torsion; supplied transverse gradient.
    s = _axis(s, 's'); n = len(s); _parameters(Ia)
    x, T, N, B = (_finite(a, 'worldline') for a in (x, T, N, B))
    Dt, Dn, Db, kappa, tau_g = (_finite(a, 'profile') for a in (Dt, Dn, Db, kappa, tau_g))
    if any(a.shape != (n, 3) for a in (x, T, N, B)) or any(a.shape != (n,) for a in (Dt, Dn, Db, kappa, tau_g)):
        raise ValueError('Worldline shape mismatch')
    if np.any(kappa <= 0):
        raise ValueError('Frenet profile requires kappa > 0')
    G = _finite(gradient_at_worldline, 'worldline gradient')
    if G.shape != (n, 3, 3):
        raise ValueError('Supply full Cartesian gradient at each sample')
    D = Dt[:, None]*T+Dn[:, None]*N+Db[:, None]*B
    displacement_components(D, T, N, B)
    dt, dn, db = (np.gradient(a, s, edge_order=2) for a in (Dt, Dn, Db))
    derivative = (dt-kappa*Dn)[:, None]*T+(dn+kappa*Dt-tau_g*Db)[:, None]*N+(db+tau_g*Dn)[:, None]*B
    return dict(D=D, derivative=derivative, omega_D=curl_displacement(G),
                tangent_residual=np.gradient(x, s, axis=0, edge_order=2)-T,
                derivative_residual=np.einsum('nij,nj->ni', G, T)-derivative)


def main():
    parser = argparse.ArgumentParser(description='Additive warp-vortex diagnostics')
    parser.add_argument('input', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--tolerance', type=float, default=1e-8)
    parser.add_argument('--identity', default='{}', help='JSON scalar parameters')
    args = parser.parse_args()
    try:
        Ia = json.loads(args.identity)
        if not isinstance(Ia, dict):
            raise ValueError('identity must be a JSON object')
        with np.load(args.input, allow_pickle=False) as data:
            result = warp_vortex_status(data['D'], (data['x'], data['y'], data['z']),
                data['kappa'], data['tau_g'], Ia,
                lambda *unused: data['sigma_residual'], lambda *unused: data['sigma_weight'],
                lambda *unused: data['envelope_level'], args.tolerance)
        # Exclusive output creation; no existing result is overwritten.
        with args.output.open('xb') as stream:
            np.savez_compressed(stream, **result)
        summary = {k: (v if np.isfinite(v) else None) for k, v in result.items() if isinstance(v, float)}
        print(json.dumps(summary, allow_nan=False))
    except (ValueError, KeyError, OSError, TypeError) as error:
        parser.error(str(error))


if __name__ == '__main__':
    main()
