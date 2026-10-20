# What background is the current ring an excitation of?

This is a reconstruction of the implemented equations, not a modification of them. Canonical real-component kinetic normalization is one half of each complex modulus squared. All quantities are dimensionless.

For the existing original field, Phi=exp(i Omega t)F with Omega=sqrt(2). The console's rest-frame equation is

    F_tt = Laplacian F - 2i Omega F_t + (1-|F|^2-4|Sigma|^2)F.

Returning to the physical Phi variable gives

    Phi_tt = Laplacian Phi - (1+|Phi|^2+4|Sigma|^2)Phi.

The physical Sigma equation is

    Sigma_tt = Laplacian Sigma + (3-4|Phi|^2-25|Sigma|^2)Sigma.

They follow from the potential

    Uphys(rho,eta) = (rho-1)^2/4 + rho + 25 eta^2/4 + (4rho-3)eta/2,
    rho=|Phi|^2 >= 0, eta=|Sigma|^2 >= 0.

The corresponding effective rotating potential is U0=Uphys-rho, exactly the effective rest functional used in the existing comparison. The full physical potential and the effective rotating potential must not be interchanged when identifying a vacuum or energy.

An exact identity is

    Uphys - 0.16 = rho^2/4 + rho/2 + 2rho eta + 25(eta-0.12)^2/4.

Every term is nonnegative for rho,eta >= 0. Therefore the zero-charge homogeneous physical ground state has Phi=0 and |Sigma|^2=0.12, with an arbitrary Sigma phase. The additive potential constant is conventional; the field values and potential differences do not depend on it.

The console instead selects Phi=exp(i sqrt(2)t), Sigma=0 at infinity. This is a charged rotating homogeneous solution. In the stated normalization it carries Phi charge density sqrt(2), kinetic energy density 1 and potential density 1. An infinite homogeneous background has infinite total energy/charge; localized excitation quantities require explicitly stated subtraction and charge matching.

This does not make its ring mathematically inconsistent. It identifies the environment: a finite-density charged background rather than the zero-charge physical ground state. Stability in that environment, if established, is not by itself derivation of an isolated vacuum particle. A claim about observed matter would need either an appropriate background-independent particle construction or a physically justified medium/background with quantitative observational consequences.

The earlier static global-Witten control uses U0 as its physical potential. It is therefore a different physical model, although the chosen effective stiffness and carrier mass match. Its comparison was explicitly not an equal-charge test of the same isolated system.

`physical_potential_audit.py` checks the force as the negative real-component potential gradient and checks the displayed nonnegative decomposition. It does not change the console, supply Authorship coupling, or prove nonlinear stability.

## A sharper binding obstruction at lambda=25

In the zero-charge ground state, the small-amplitude Phi mass squared is m_Phi^2=1+4(0.12)=1.48. A second exact decomposition is

    Uphys - Umin - 0.74 rho
        = 25/4 (eta - 0.12 + 4rho/25)^2 + 0.09 rho^2 >= 0.

For finite-energy fields over that ground state, the canonical energy therefore obeys

    E_relative >= integral [(|Phi_t|^2 + m_Phi^2 |Phi|^2)/2]
               >= m_Phi |Q_Phi|,
    Q_Phi = integral Im(Phi* Phi_t).

The first inequality discards nonnegative spatial gradients, Sigma kinetic energy and the displayed potential squares. The second follows pointwise from the arithmetic/geometric mean inequality and |Im(Phi* Phi_t)| <= |Phi| |Phi_t|, then the triangle inequality on the integrated charge. A dilute Phi charge distribution approaches the same free threshold. Consequently the current lambda=25 potential supplies no below-free-threshold globally bound Phi-charge soliton in this vacuum sector. This is not a theorem excluding every excited/metastable structure, does not exclude medium-supported rings, and is not a comparison at nonzero additional Sigma charge.

## Existing lambda=12 family as a separate question

The code archive already explores lambda=12 in a different rotating-background search. That is not a vacuum solution. However the same physical potential family can be audited without changing any console equations. For general lambda>0 the true vacuum has rho=0, eta=3/lambda, Umin=1/4-9/(4lambda), and m_Phi^2=1+12/lambda.

For rho>0, minimizing 2(Uphys-Umin)/rho over eta>=0 gives two ranges: eta=(3-4rho)/lambda for rho<=3/4, and eta=0 for rho>=3/4. Minimizing subsequently over rho yields

    mu_star^2 = 1+3/sqrt(lambda)       for lambda<16,
                1+12/lambda          for lambda>=16.

For lambda>=16 the latter is an infimum at rho approaching zero. For lambda<16 the minimum occurs at rho=3/sqrt(lambda), eta=0. At lambda=12, mu_star^2=1.8660254 < m_Phi^2=2. This permits a pointwise subthreshold well; gradients, fixed-charge solution construction and stability remain unresolved. It is a mathematically motivated next soliton question, not a demonstrated particle or license to adjust coefficients to force agreement with observed matter.
