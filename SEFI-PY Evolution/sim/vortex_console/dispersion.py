"""Exact rest-frame far-field response; does not modify vortex evolution."""
import math
import numpy as np

MODEL_VERSION = 'far-field-dispersion-v1'

def frequencies(k, omega):
    k = np.asarray(k, dtype=float)
    if not np.isfinite(k).all() or np.any(k < 0):
        raise ValueError('Wave numbers must be finite and nonnegative')
    if not math.isfinite(omega) or omega < 0:
        raise ValueError('Background frequency must be finite and nonnegative')
    a = 1 + 2 * omega**2
    root = np.sqrt(a*a + 4*omega**2*k*k)
    # Rationalization preserves precision on the gapless branch near k=0.
    low2 = k*k*(k*k+2)/(k*k+a+root)
    return np.sqrt(low2), np.sqrt(k*k+a+root)

def linear_matrix(k, omega):
    return np.array([[0,0,1,0], [0,0,0,1],
                     [-(k*k+2),0,0,2*omega], [0,-k*k,-2*omega,0]], dtype=float)

def comparison():
    k = np.linspace(0, 3, 241)
    models = []
    for name, omega in [('Static global-Witten control', 0.),
                        ('Rotating console background', math.sqrt(2))]:
        low, high = frequencies(k, omega)
        speed = 1/math.sqrt(1+2*omega**2)
        phase_speed = np.divide(low, k, out=np.full_like(k, speed), where=k != 0)
        models.append(dict(name=name, omega=omega, low=low.tolist(), high=high.tolist(),
                           phase_speed=phase_speed.tolist(), long_wave_speed=speed))
    return dict(model_version=MODEL_VERSION, k=k.tolist(), models=models,
                speed_reduction_percent=100*(1-1/math.sqrt(5)),
                scope='Linear homogeneous far-field response at rest; not nonlinear vortex evolution.',
                matching='Same effective amplitude stiffness and carrier mass; different physical potentials and background charges.',
                interpretation='Not evidence of unique SEFI physics, improved stability, or a changed causal light speed.',
                reference='https://arxiv.org/abs/1909.01950')
