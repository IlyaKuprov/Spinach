# kernel/utilities/dirdiff.m

## Purpose

Computes directional derivatives of the matrix exponential, implementing Equation 11 of Najfeld and Havel and Equation 16 of Goodwin and Kuprov. The function returns the propagator and its derivatives with respect to perturbations of the Hamiltonian along one or more specified directions.

Source: [kernel/utilities/dirdiff.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/dirdiff.m)

## Behaviour

The function builds an auxiliary block matrix of dimension `N` by `N` blocks, each block having the size of `A`. All diagonal blocks are set to `A`. The superdiagonal blocks are set to the differentiation direction or directions: if `B` is a cell array, block `{n,n+1}` receives `B{n}` for `n=1..N-1`; if `B` is a single matrix, every superdiagonal block receives `B`. All other blocks are initialised as sparse zero matrices of the same size as `A`.

Before exponentiation, the propagator tolerance is tightened by setting `spin_system.tols.prop_chop` to `1e-14`. The auxiliary matrix is converted with `cell2mat` and exponentiated over time `T` using the `propagator` function, yielding `exp(-1i*A*T)` in the top-left block together with the derivative blocks.

Directional derivatives are extracted from the first block row: output element `n` is `factorial(n-1)` times the block `auxmat(1:size(A,1), (1:size(A,2))+size(A,2)*(n-1))`. The result is the cell array `{D0,D1,D2,...}` of Equation 18 in Goodwin and Kuprov, where `D0` is the propagator itself and subsequent entries are its directional derivatives. Using `N=2` yields the propagator and its first derivative.

Input consistency is enforced by an internal `grumble` subroutine, which raises errors when:

- `N` is not a real scalar integer greater than 1;
- `B` is a cell array whose number of elements does not equal `N-1`;
- `A` or any `B` matrix is not a numeric square matrix;
- `T` is not a real numeric scalar.

## Inputs and outputs

**Inputs**

- `spin_system` — spin system object supplying tolerances used by the propagator call.
- `A` — Hamiltonian at the reference point, corresponding to the `exp(-1i*A*T)` propagator; must be a numeric square matrix.
- `B` — differentiation direction (a single square matrix) or directions (a cell array of square matrices); when a cell array, the number of matrices must equal `N-1`.
- `T` — time used in `exp(-1i*A*T)`; must be a real scalar.
- `N` — block dimension of the auxiliary matrix; must be a real integer greater than 1. Use `N=2` to obtain the propagator and its first derivative.

**Outputs**

- `D` — cell array of matrices `{D0,D1,D2,...}` corresponding to Equation 18 in Goodwin and Kuprov.

## References

- Najfeld, I.; Havel, T. F. Derivatives of the matrix exponential and their computation. *Advances in Applied Mathematics*. [https://doi.org/10.1006/aama.1995.1017](https://doi.org/10.1006/aama.1995.1017)
- Goodwin, D. L.; Kuprov, I. Auxiliary matrix formalism for interaction picture transformations of pulsed and continuous wave calculations. *The Journal of Chemical Physics*. [https://doi.org/10.1063/1.4928978](https://doi.org/10.1063/1.4928978)
- Spinach documentation: [dirdiff.m](https://spindynamics.org/wiki/index.php?title=dirdiff.m)
