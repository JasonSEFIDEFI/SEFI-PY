# CP098-T90C-R29AR — Universal positivity proof

## Domain

x > 0
u = y - 31/32 > 0
q = h^2 > 0

## Established prerequisites

R29AG:
The exact reduced resultant is R = x^4 G,
where G is cubic in q.

R29AH and R29AO:
B0, B2, B3 are strictly positive throughout x,u > 0.

R29AL:
B1 has exactly one positive x-root r(u).
B1 >= 0 for 0 < x <= r(u), and B1 < 0 afterward.

R29AQ:
The cubic discriminant Delta has exactly one
positive x-root d(u).

R29AS:
The exact discriminant has a strictly positive
constant coefficient and strictly negative
leading coefficient for every u > 0.
Combined with R29AQ, Delta is positive before
d(u) and negative afterward.

R29AO:
Delta < 0 wherever B1 = 0.

## Root ordering

For every u > 0, Delta(r(u),u) < 0.
The unique positive discriminant root therefore
satisfies 0 < d(u) < r(u).

## Cubic positivity

G(q) = B0 + B1*q + B2*q^2 + B3*q^3.

If 0 < x <= r(u), all terms are nonnegative
and B0 > 0, hence G(q) > 0 for q > 0.

If x > r(u), then Delta < 0.
The cubic has exactly one distinct real root.
Because B3 > 0 and B0 > 0, that root is
strictly negative. Hence G(q) > 0 for q > 0.

Therefore G(x,u,q) > 0 throughout x,u,q > 0.

## Resultant consequence

R = x^4 G > 0 throughout the stated domain.
Thus the reduced polynomials PD and PT have
no common complex root in b at any admissible
parameter point in this domain.

## Scope and stop-gates

This establishes positivity for the exact
reduced polynomial system and stated domain.

It does not independently establish equivalence
to every unreduced original constraint, cover
unexamined boundary branches, or prove full
Dirac closure.

Physical validation is not claimed.
