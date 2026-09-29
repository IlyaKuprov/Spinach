# kernel/utilities/expmint2.m

## Purpose

`expmint2` computes the nested matrix exponential double integral

`Integrate[expm(-i*A*(T-t))*B*Integrate[expm(-i*C*(t-x))*D*expm(-i*E*x),{x,0,t}],{t,0,T}]`

as documented in the function header. The result corresponds to the (1,3) block of the exponential of an auxiliary block matrix, following the method of Van Loan (http://dx.doi.org/10.1109/TAC.1978.1101743).

## Behaviour

- Syntax: `I=expmint2(spin_system,A,B,C,D,E,T)`.
- The function first runs an internal consistency check (`grumble`) on all arguments.
- A zero filler block `Z` is created as a sparse matrix with the dimensions of `A`.
- The auxiliary matrix is assembled as a 3-by-3 block matrix:

  `auxmat = [A  -1i*B,     Z; Z      C  -1i*D; Z      Z      E]`

- The auxiliary matrix is exponentiated over the interval `T` using `propagator(spin_system,auxmat,T)`.
- Block extractors `BE1=[speye(size(A))  Z  Z]` and `BE3=[Z; Z; speye(size(A))]` are built, and the integral is extracted as `I=BE1*P*BE3`, i.e. the (1,3) block of the propagated auxiliary matrix.
- Consistency enforcement (`grumble`) requires:
  - All of `A`, `B`, `C`, `D`, `E`, `T` to be numeric, otherwise `'all arguments must be numeric.'` is raised.
  - `A`, `B`, `C`, `D`, `E` to be matrices, otherwise `'A, B, C, D, E must be matrices.'` is raised.
  - All matrices to be square, otherwise `'all matrices must be square.'` is raised.
  - All matrices to have the same dimension, otherwise `'all matrices must have the same dimension.'` is raised.
  - `T` to be a real scalar, otherwise `'T must be a real scalar.'` is raised.

## Inputs and outputs

Inputs:

- `spin_system` — spin system object passed through to `propagator`.
- `A`, `B`, `C`, `D`, `E` — square matrices of the same dimension.
- `T` — upper limit of the outer integral; a real scalar.

Output:

- `I` — the nested matrix exponential double integral as defined above.

## References

- Source: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/expmint2.m
- Wiki: https://spindynamics.org/wiki/index.php?title=expmint2.m
- C. F. Van Loan, computing integrals involving the matrix exponential, IEEE Transactions on Automatic Control, http://dx.doi.org/10.1109/TAC.1978.1101743
