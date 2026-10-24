import numpy as np
from audit import exact_square_identity, potential, minimized_ratio, analytic_audit, solve_profile


def test_exact_identity_and_lambda_transition():
    assert all(exact_square_identity(lam) for lam in (1, 4, 12, 16, 25, 100))
    for lam in (4., 12., 15., 16., 25.):
        rho = np.geomspace(1e-6, 30, 1000)
        eta = np.geomspace(1e-7, 20, 1000)
        lhs = potential(rho, eta, lam)-(1+12/lam)*rho/2
        rhs = lam/4*(eta-3/lam+4*rho/lam)**2+(.25-4/lam)*rho**2
        np.testing.assert_allclose(lhs, rhs, atol=1e-10)
        if lam >= 16:
            assert np.min(lhs) >= -1e-10
        else:
            x = 3/np.sqrt(lam)
            assert minimized_ratio(x, lam) < 1+12/lam


def test_vacuum_hessian_independent_real_components():
    lam = 12.; v = np.sqrt(3/lam)
    vacuum = np.array([0., 0., v, 0.])
    def U(x):
        return potential(x[0]**2+x[1]**2, x[2]**2+x[3]**2, lam)
    epsilon = 1e-4
    hessian = np.zeros((4, 4))
    for i in range(4):
        for j in range(4):
            ei = np.eye(4)[i]*epsilon; ej = np.eye(4)[j]*epsilon
            hessian[i, j] = (U(vacuum+ei+ej)-U(vacuum+ei-ej)-U(vacuum-ei+ej)+U(vacuum-ei-ej))/(4*epsilon**2)
    np.testing.assert_allclose(hessian, np.diag([2, 2, 6, 0]), atol=1e-6)


def test_numerical_potential_minimization():
    for row in analytic_audit()['coupling_audit']:
        assert abs(row['numeric_infimum']-row['infimum_ratio']) < 1e-7


def test_fresh_profile_equations_and_charge():
    result, model, f, s = solve_profile(100000., .4, tolerance=1e-7)
    assert result['usable_stationarity_1e-7']
    # Independently assemble finite-volume gradient, not model.gradient.
    for field, source in ((f, 0.), (s, model.vac)):
        flux = np.zeros(len(field)+1)
        faces = np.arange(len(field)+1)*model.h
        flux[1:-1] = 4*np.pi*faces[1:-1]**2*np.diff(field)/model.h
        flux[-1] = 4*np.pi*model.L**2*(source-field[-1])/(model.h/2)
        laplace = np.diff(flux)/model.V
        force = ((result['omega']**2-1-f*f-4*s*s)*f if source == 0 else
                 (3-4*f*f-model.lam*s*s)*s)
        assert np.max(abs(laplace+force)) < 1e-7
    assert abs(result['omega']*np.sum(model.V*f*f)-100000.) < 1e-8


def test_dynamic_linearization_independent_rhs():
    from dynamics import matrices
    from audit import Radial
    model = Radial(Q=30., h=1., L=10.)
    rng = np.random.default_rng(172)
    f = rng.uniform(.1, .7, len(model.r)); s = rng.uniform(.1, .5, len(model.r))
    A, _, _, _ = matrices(model, f, s)
    omega = model.Q/np.sum(model.V*f*f)
    root = np.sqrt(model.V)
    perturbation = rng.normal(size=6*len(f))
    u, v, a, pu, pv, pa = np.split(perturbation, 6)
    def rhs(epsilon):
        F = f+epsilon*(u+1j*v)/root
        S = s+epsilon*a/root
        P = epsilon*(pu+1j*pv)/root
        sigma_velocity = epsilon*pa/root
        acc_phi = -model.K(F)/model.V-2j*omega*P+(omega**2-1-abs(F)**2-4*abs(S)**2)*F
        acc_sigma = -(model.K(S)-model.b)/model.V+(3-4*abs(F)**2-model.lam*abs(S)**2)*S
        return np.r_[root*P.real, root*P.imag, root*sigma_velocity,
                     root*acc_phi.real, root*acc_phi.imag, root*acc_sigma.real]
    finite = (rhs(1e-5)-rhs(-1e-5))/(2e-5)
    np.testing.assert_allclose(finite, A@perturbation, rtol=1e-7, atol=1e-7)


def test_linearized_charge_and_phase_symmetry():
    from dynamics import matrices
    result, model, f, s = solve_profile(100000., .8)
    assert result['usable_stationarity_1e-7']
    A, local, constrained, _ = matrices(model, f, s)
    n = len(f); root_f = np.sqrt(model.V)*f
    phase = np.zeros(6*n); phase[n:2*n] = root_f
    charge = np.zeros(6*n); charge[:n] = 2*result['omega']*root_f; charge[4*n:5*n] = root_f
    gf, _, _ = model.gradient(f, s)
    phase_residual = np.zeros(6*n); phase_residual[4*n:5*n] = np.sqrt(model.V)*gf
    charge_residual = np.zeros(6*n); charge_residual[n:2*n] = np.sqrt(model.V)*gf
    np.testing.assert_allclose(A@phase, phase_residual, atol=1e-11)
    np.testing.assert_allclose(charge@A, charge_residual, atol=1e-11)
    assert np.linalg.norm(A@phase)/np.linalg.norm(phase) < 1e-7
    assert np.linalg.norm(charge@A)/np.linalg.norm(charge) < 1e-7
    # The elimination term is positive rank one, not an ad hoc stabilizer.
    update = constrained-local
    assert np.linalg.eigvalsh(update).min() > -1e-10
    assert np.linalg.matrix_rank(update, tol=1e-9) == 1
