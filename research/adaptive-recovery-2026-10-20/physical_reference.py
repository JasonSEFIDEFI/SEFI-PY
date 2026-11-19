"""Physical oscillator reference and unconditional metrics. Jason D. Dutton.

Units: [a,a^dagger]=1, beta=x+iy, q=(a+a^dagger)/sqrt(2).
The payload code is span{|0>,|1>}; it is an encoding, not a protecting QEC code.
Analytic overlaps refer to the infinite oscillator, avoiding cutoff bias in scores.
"""
import math
import numpy as np


def annihilation(cutoff):
    if not isinstance(cutoff, int) or cutoff < 2:
        raise ValueError("cutoff must be an integer >= 2")
    return np.diag(np.sqrt(np.arange(1, cutoff)), 1).astype(complex)


def coherent_state(alpha, cutoff):
    """Normalized projected coherent state, with retained norm reported separately."""
    if not np.isfinite(alpha):
        raise ValueError("alpha must be finite")
    annihilation(cutoff)
    amplitudes = np.empty(cutoff, complex)
    amplitudes[0] = np.exp(-abs(alpha)**2/2)
    for n in range(1, cutoff):
        amplitudes[n] = amplitudes[n-1]*alpha/math.sqrt(n)
    retained = float(np.vdot(amplitudes, amplitudes).real)
    if retained <= 0:
        raise ValueError("coherent state underflow")
    return amplitudes/math.sqrt(retained), retained


def displacement(beta, cutoff):
    """Unitary finite-oscillator approximation exp(beta a^dagger-beta* a)."""
    if not np.isfinite(beta):
        raise ValueError("beta must be finite")
    a = annihilation(cutoff)
    hermitian = 1j*(beta*a.conj().T-beta.conjugate()*a)
    values, vectors = np.linalg.eigh(hermitian)
    return (vectors*np.exp(-1j*values))@vectors.conj().T


def code_overlap(beta):
    """Exact V^dagger D(beta) V, V=(|0>,|1>)."""
    beta = np.asarray(beta, complex)
    r2 = np.abs(beta)**2
    e = np.exp(-r2/2)
    m = np.empty(beta.shape+(2, 2), complex)
    m[..., 0, 0] = e
    m[..., 0, 1] = -beta.conjugate()*e
    m[..., 1, 0] = beta*e
    m[..., 1, 1] = (1-r2)*e
    return m


def metrics(residual, tolerance=0.25):
    """Entanglement and six-state mean fidelity; all leaked mass counts as failure."""
    residual = np.asarray(residual, float)
    if residual.ndim != 2 or residual.shape[1] != 2 or not np.isfinite(residual).all():
        raise ValueError("residual must be a finite n by 2 matrix")
    if tolerance <= 0:
        raise ValueError("tolerance must be positive")
    r2 = np.sum(residual**2, axis=1)
    e = np.exp(-r2)
    fe = e*(1-r2/2)**2
    survival = e*(1+r2*r2/2)
    # Haar qubit average: (Tr(M^dagger M)+|Tr M|^2)/(2*3).
    # This equals the exact average over the six Pauli eigenstates.
    state_fidelity = (survival+2*fe)/3
    return {"entanglement_fidelity": fe, "state_fidelity": state_fidelity,
            "survival": survival, "accuracy": (r2 <= tolerance**2).astype(float),
            "residual_mse": r2}


def erasure_kraus(beta):
    """CPTP qubit-to-qutrit decoded map: code block plus an orthogonal erasure flag."""
    m = code_overlap(beta)
    if m.shape != (2, 2):
        raise ValueError("single beta required")
    k = np.zeros((3, 2), complex); k[:2] = m
    complement = np.eye(2)-m.conj().T@m
    values, vectors = np.linalg.eigh(complement)
    if values.min() < -1e-12:
        raise ValueError("nonphysical code block")
    result = [k]
    for value, vec in zip(values, vectors.T):
        leak = np.zeros((3, 2), complex)
        leak[2] = math.sqrt(max(0, value))*vec.conjugate()
        result.append(leak)
    return result


def choi(kraus):
    """Unnormalized Choi matrix, column-major vectorization; trace equals input dim."""
    vecs = [k.reshape(-1, order="F") for k in kraus]
    return sum(np.outer(v, v.conjugate()) for v in vecs)
