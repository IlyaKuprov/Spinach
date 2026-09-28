# kernel/utilities/expmint.m

- Signature: `R=expmint(spin_system,A,B,C,T)`

## Purpose

Computes matrix exponential integrals of the following general type: Integrate[expm(-i*A*t)*B*expm(i*C*t),{t,0,T}] Matrix A must be Hermitian. For further info see the paper by Char- les van Loan (http://dx.doi.org/10.1109/TAC.1978.1101743). Syntax: R=expmint(spin_system,A,B,C,T)

## Physical / mathematical content

- Evaluates the matrix integral of `exp(-1i*A*t)*B*exp(1i*C*t)` over `0 <= t <= T`, with Hermitian `A`.

## Numerical / algorithmic content

- For nonzero `T` and nonzero `B`, forms an auxiliary block matrix, computes its matrix exponential, and extracts the integral. Returns a sparse zero matrix when `T` is zero or `B` has no nonzero entries.

## Parameters / inputs

- A,B,C -the three matrices involved in the integral
- T -integration time
- Output:
- R -the resulting integral
- Note: the auxiliary matrix method is massively faster than either
- commutator series or diagonalisation.
- Note: this is the most memory-intensive stage in a lot of calcula-
- tions; memory recycling is aggressive.

## Implementation structure

- Computes matrix exponential integrals of the following general type:
- Integrate[expm(-i*A*t)*B*expm(i*C*t),{t,0,T}]
- Matrix A must be Hermitian. For further info see the paper by Char-
- les van Loan (http://dx.doi.org/10.1109/TAC.1978.1101743). Syntax:
- R=expmint(spin_system,A,B,C,T)
- A,B,C -the three matrices involved in the integral
- T -integration time
- Output:
- R -the resulting integral
- Note: the auxiliary matrix method is massively faster than either
- commutator series or diagonalisation.
- Note: this is the most memory-intensive stage in a lot of calcula-
