# kernel/utilities/expmint2.m

- Signature: `I=expmint2(spin_system,A,B,C,D,E,T)`

## Purpose

Computes the nested matrix exponential double integral: Integrate[expm(-i*A*(T-t))*B* Integrate[expm(-i*C*(t-x))*D*expm(-i*E*x),{x,0,t}],{t,0,T}] This corresponds to the (1,3) block of the exponential of the auxiliary matrix (http://dx.doi.org/10.1109/TAC.1978.1101743). Syntax: I=expmint2(spin_system,A,B,C,D,E,T)

## Physical / mathematical content

- Computes the nested matrix-exponential double integral stated above, involving `A`, `B`, `C`, `D`, `E`, and the interval from zero to `T`.

## Numerical / algorithmic content

- Builds a 3-by-3 upper block-triangular auxiliary matrix with diagonal blocks `A`, `C`, and `E`, calls `propagator` for time `T`, and extracts the `(1,3)` block.

## Parameters / inputs

- A,B,C,D,E -square matrices
- T -upper limit of the outer integral
- Output:
- I -the integral as above

## Implementation structure

- Forms the auxiliary matrix `[A, -1i*B, 0; 0, C, -1i*D; 0, 0, E]`, evaluates its propagator, and extracts the `(1,3)` block using left and right block-selection matrices.
