"""Reproduce selected Time manuscript results. Author: Jason D. Dutton."""
import math
from fractions import Fraction
from core.time_dynamics import (
    IdentityMetric, IdentityModel, Coupling, exact_linear_ratios,
    recurrence_cycles, floquet_monodromy, floquet_diagnostics,
)


def main():
    print('SEFI-DEFI-GWFM: dimensionless classical model, not empirical validation')
    ratios=exact_linear_ratios(Fraction(24,25))
    print('r=24/25: exact frequency ratios:',ratios)
    print('Full-state recurrence:',recurrence_cycles(ratios),'Origin cycles (10*pi)')
    nonlinear=IdentityModel(IdentityMetric(),Coupling.quadratic(.5),Coupling.linear(0))
    print('Nonlinear barrier:',nonlinear.effective_potential(math.sqrt(2),1,0))
    for omega in (math.sqrt(2),2.2):
        matrix=floquet_monodromy(omega)
        result=floquet_diagnostics(matrix,math.tau/omega)
        print('Drive:',omega,'monodromy:',matrix)
        print('Result:',result)


if __name__=='__main__':
    main()
