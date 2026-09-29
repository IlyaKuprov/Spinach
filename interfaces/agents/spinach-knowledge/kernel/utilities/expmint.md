# kernel/utilities/expmint.m

## Purpose

Computes matrix exponential integrals of the general type

`Integrate[expm(-i*A*t)*B*expm(i*C*t),{t,0,T}]`

using the auxiliary matrix method of Charles van Loan. The matrix `A` must be Hermitian.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/expmint.m>

## Behaviour

- Syntax: `R=expmint(spin_system,A,B,C,T)`.
- Consistency is enforced first: `A`, `B` and `C` must be numeric matrices of the same dimension, `A` must be Hermitian, and `T` must be a real scalar; violations raise errors.
- Zero-argument shortcut: if `T==0` or `B` has no nonzero entries, the function returns an empty sparse matrix of the size of `B` via `spalloc`.
- Otherwise the auxiliary matrix `[-A, 1i*B; 0*A, -C]` is built, block extraction multipliers `BE1` and `BE2` are formed with `spdiags`, and `A`, `B`, `C` are cleared from memory.
- The auxiliary matrix is exponentiated over time `T` by calling `propagator(spin_system,auxmat,T)`.
- The user is informed via `report(spin_system,'processing auxiliary matrix blocks...')`.
- Blocks are extracted as `P=(BE1'*auxmat*BE1)'` and `Q=(BE1'*auxmat*BE2)`, after which the auxiliary matrix and extraction multipliers are cleared.
- The result is `R=clean_up(spin_system,P*Q,spin_system.tols.prop_chop)`, after which `P` and `Q` are cleared.
- The header notes that the auxiliary matrix method is massively faster than either commutator series or diagonalisation, and that this is the most memory-intensive stage in many calculations, with aggressive memory recycling.

## Inputs and outputs

Inputs:

- `spin_system` — Spinach system object passed through to `propagator`, `report` and `clean_up`.
- `A`, `B`, `C` — the three matrices involved in the integral; they must be numeric matrices of identical dimension, and `A` must be Hermitian.
- `T` — integration time; a real scalar.

Output:

- `R` — the resulting integral.

## References

- C. van Loan, paper on computing integrals involving matrix exponentials: <http://dx.doi.org/10.1109/TAC.1978.1101743>
- Spinach Wiki page: <https://spindynamics.org/wiki/index.php?title=expmint.m>
