"""Internal symmetry/torque audit on prescribed worldline geometry.

This does not construct bulk fields or solve the worldline geometry equation.
"""
import json
from pathlib import Path
from fractions import Fraction
import numpy as np


def symmetry_conditions(h, k):
    h, k = list(map(Fraction, h)), list(map(Fraction, k))
    if len(h) != 4 or len(k) != 4 or min(h+k) <= 0:
        raise ValueError('four positive kinetic weights and stiffnesses required')
    return {'ordinary_origin_centered_pair_rotations': [h[i] == h[j] and k[i] == k[j] for i, j in ((0, 1), (2, 3))],
            'kinetically_normalized_pair_rotations': [k[i]/h[i] == k[j]/h[j] for i, j in ((0, 1), (2, 3))],
            'squared_frequencies': [str(ki/hi) for hi, ki in zip(h, k)]}


def potential(q, h, k, curvature, torsion, coefficients):
    I = np.asarray(q)/np.sqrt(h)
    c1, d1, c2, d2 = coefficients
    return float(.5*np.dot(np.asarray(k)/h, q*q)+(c1+d1*I[1])*curvature**2*I[2]
                 +(c2+d2*I[1])*torsion**2*I[3])


def acceleration(q, h, k, curvature, torsion, coefficients):
    I = np.asarray(q)/np.sqrt(h)
    c1, d1, c2, d2 = coefficients
    force = np.array([0., (d1*curvature**2*I[2]+d2*torsion**2*I[3])/np.sqrt(h[1]),
                      (c1+d1*I[1])*curvature**2/np.sqrt(h[2]),
                      (c2+d2*I[1])*torsion**2/np.sqrt(h[3])])
    return -np.asarray(k)/h*q-force


def charges(q, velocity):
    return np.array([q[0]*velocity[1]-q[1]*velocity[0], q[2]*velocity[3]-q[3]*velocity[2]])


def torque(q, h, k, curvature, torsion, coefficients):
    I = q/np.sqrt(h); omega2 = np.asarray(k)/h
    c1, d1, c2, d2 = coefficients
    return np.array([(omega2[0]-omega2[1])*q[0]*q[1]
                     -q[0]/np.sqrt(h[1])*(d1*curvature**2*I[2]+d2*torsion**2*I[3]),
                     (omega2[2]-omega2[3])*q[2]*q[3]
                     +q[3]*curvature**2*(c1+d1*I[1])/np.sqrt(h[2])
                     -q[2]*torsion**2*(c2+d2*I[1])/np.sqrt(h[3])])


def findings():
    h = np.array([1., 4., 9., 16.]); k = np.array([2., 8., 27., 48.])
    q = np.array([1., 2., 3., 4.])
    cases = []
    for name, coefficients in [('uncoupled', [0., 0., 0., 0.]),
                               ('constant_geometry_coupling', [1., 0., 2., 0.]),
                               ('authorship_dependent_coupling', [1., .5, 2., -.25])]:
        cases.append({'name': name, 'coefficients_c1_d1_c2_d2': coefficients,
                      'dJ_ds': torque(q, h, k, 2., 3., coefficients).tolist()})
    return {'scope': 'Prescribed kappa=2, torsion=3, q=(1,2,3,4); exact symmetry conditions and illustrative torques',
            'h': h.tolist(), 'k': k.tolist(), 'conditions': symmetry_conditions([1, 4, 9, 16], [2, 8, 27, 48]),
            'cases': cases,
            'interpretation': 'Constant couplings preserve OA rotations but generally torque origin-centered SW rotations. Affine shifted symmetries may survive on fixed constant geometry.'}


if __name__ == '__main__':
    path = Path(__file__).with_name('internal-symmetry-results.json')
    if path.exists():
        raise FileExistsError(path)
    path.write_text(json.dumps(findings(), indent=2)+'\n')
