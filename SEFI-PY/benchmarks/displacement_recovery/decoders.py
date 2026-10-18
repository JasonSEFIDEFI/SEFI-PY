"""Measurement-only estimators, matched information. Jason D. Dutton.

Geometry = positive quadratic phase-space metric learned from calibration.
It is ordinary Gaussian regularization, motivated by reference-geometry distances.
No claim of a new inference algorithm or a unique SEFI deduction is made.
"""
import numpy as np


def covariance_from_calibration(calibration):
    x = np.asarray(calibration, float)
    if x.ndim != 2 or x.shape[1] != 2 or len(x) < 3 or not np.isfinite(x).all():
        raise ValueError("at least three finite 2D calibration injections required")
    return np.cov(x, rowvar=False, ddof=1)


class MaximumLikelihood:
    name = "ML"
    def predict(self, observations):
        return np.asarray(observations, float).copy()


class GeometryDecoder:
    """min_z (y-z)^T R^-1(y-z)+z^T Sigma^-1 z, with zero prior mean.

    Zero covariance is the well-defined zero-displacement prior limit.
    Model setup accepts calibration covariance and measurement noise only.
    """
    name = "Geometry"
    def __init__(self, sigma, measurement_covariance):
        s = np.asarray(sigma, float); r = np.asarray(measurement_covariance, float)
        if s.shape != (2, 2) or r.shape != (2, 2) or not np.isfinite(s+r).all():
            raise ValueError("finite 2 by 2 covariance matrices required")
        if not np.allclose(s,s.T) or not np.allclose(r,r.T):
            raise ValueError("covariances must be symmetric")
        if np.linalg.eigvalsh(s).min() < -1e-12 or np.linalg.eigvalsh(r).min() <= 0:
            raise ValueError("prior PSD and measurement PD required")
        self.weight = np.linalg.solve(s+r, s).T
    def predict(self, observations):
        return np.asarray(observations, float)@self.weight.T


class ConventionalPosterior:
    """Independent precision-form Gaussian posterior mean, same fitted prior."""
    name = "Conventional posterior"
    def __init__(self, sigma, measurement_covariance):
        if np.linalg.norm(sigma) == 0:
            self.weight = np.zeros((2, 2))
        else:
            precision = np.linalg.inv(measurement_covariance)
            self.weight = np.linalg.solve(np.linalg.inv(sigma)+precision, precision)
    def predict(self, observations):
        return np.asarray(observations, float)@self.weight.T


class AffineCoordinates:
    name = "Affine coordinates"
    def __init__(self, sigma, r, transform):
        self.transform = np.asarray(transform, float)
        self.inverse = np.linalg.inv(self.transform)
        self.decoder = GeometryDecoder(self.transform@sigma@self.transform.T,
                                       self.transform@r@self.transform.T)
    def predict(self, observations):
        transformed = np.asarray(observations)@self.transform.T
        return self.decoder.predict(transformed)@self.inverse.T


class NonlinearCoordinates:
    """Bijective sinh/arcsinh chart round trip, with physical metric pulled back.

    This is an information-preserving representation check, not evidence of
    intrinsic nonlinear-geometric robustness or invariant coordinate MAP density.
    """
    name = "Nonlinear coordinates"
    def __init__(self, decoder):
        self.decoder = decoder
    def predict(self, observations):
        encoded = np.arcsinh(np.asarray(observations))
        return self.decoder.predict(np.sinh(encoded))


class ShuffledFeatures:
    name = "Shuffled features"
    def __init__(self, decoder, permutation):
        self.decoder = decoder; self.permutation = permutation
    def predict(self, observations):
        return self.decoder.predict(np.asarray(observations)[self.permutation])


def make_decoders(sigma, r, transform, rng, trials):
    geometry = GeometryDecoder(sigma, r)
    isotropic = np.eye(2)*np.trace(sigma)/2
    angle = rng.uniform(0, 2*np.pi)
    rot = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    wrong_sigma = rot@sigma@rot.T
    return {
        "No recovery": GeometryDecoder(np.zeros((2, 2)), r),
        "ML": MaximumLikelihood(),
        "Geometry": geometry,
        "Conventional posterior": ConventionalPosterior(sigma, r),
        "No prior distance": MaximumLikelihood(),
        "No anisotropy": GeometryDecoder(isotropic, r),
        "Random metric orientation": GeometryDecoder(wrong_sigma, r),
        "Shuffled features": ShuffledFeatures(geometry, rng.permutation(trials)),
        "Affine coordinates": AffineCoordinates(sigma, r, transform),
        "Nonlinear coordinates": NonlinearCoordinates(geometry)
    }
