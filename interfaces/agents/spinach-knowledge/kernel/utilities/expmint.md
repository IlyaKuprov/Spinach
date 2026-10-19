# kernel/utilities/expmint.m

## Purpose

Computes matrix exponential integrals of the general type

`Integrate[expm(-i*A*t)*B*expm(i*C*t),{t,0,T}]`

using the auxiliary matrix method of Charles van Loan. The matrix `A` must be Hermitian.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/expmint.m>

## Behaviour

Use `R=expmint(spin_system,A,B,C,T)` for the integral above. The matrices must be numeric and have identical dimensions, `A` must be Hermitian, and `T` must be a real scalar; invalid arguments raise errors. Zero integration time or a zero `B` gives a zero sparse result with the dimensions of `B`.

The auxiliary-matrix method provides an alternative to commutator expansions or diagonalisation. It doubles the matrix dimension and can be the most memory-intensive part of a calculation, so the propagation tolerances and sparse/full storage policy in `spin_system` matter. Intermediate storage is recycled aggressively; retain sparse inputs when appropriate rather than assuming that a sparse input guarantees a small temporary product.

## Inputs and outputs

Inputs:

- `spin_system` — Spinach configuration controlling propagation tolerances, cleanup and reporting.
- `A`, `B`, `C` — the three matrices involved in the integral; they must be numeric matrices of identical dimension, and `A` must be Hermitian.
- `T` — integration time; a real scalar.

Output:

- `R` — the resulting integral.

## Temporary memory

Scaling-and-squaring and the final back-rotation use `multiclean` to discard insignificant product entries panel by panel rather than retaining the complete unpruned product. The existing cleanup tolerance and whole-matrix storage policy remain in effect. This reduces temporary fill-in at a possible runtime cost; it does not reduce the auxiliary dimension or change the integral formula.

## References

- C. van Loan, paper on computing integrals involving matrix exponentials: <http://dx.doi.org/10.1109/TAC.1978.1101743>
- Spinach Wiki page: <https://spindynamics.org/wiki/index.php?title=expmint.m>
