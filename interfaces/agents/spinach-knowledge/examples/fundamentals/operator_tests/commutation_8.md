# examples/fundamentals/operator_tests/commutation_8.m

- Signature: `commutation_8()`

## Purpose

Checks Hilbert-to-Liouville operator-action identities and the mapping between first-rank Stevens operators and Pauli Cartesian operators.

## Mathematical content

For a matrix `R`, left multiplication by `H`, right multiplication by `H`, their commutator action, and their anticommutator action have distinct Liouville-space representations. This function compares those representations with the corresponding explicit Hilbert-space products after MATLAB column-major vectorisation. It separately checks the first-rank Stevens-to-Pauli correspondence at several multiplicities.

## Callable context and model

Call the zero-input MATLAB function `commutation_8()` from a Spinach checkout with the project functions on the MATLAB path. It returns no values, prints a success message for each test group, and raises an error when a comparison exceeds the source threshold. The Hilbert-Liouville test uses independent unseeded complex random 6-by-6 matrices `H` and `R`; they are not explicitly symmetrised or constrained to represent Hermitian physical observables.

## Checks encoded in the source

- Forms `L_left=hilb2liouv(H,'left')`, `L_right=hilb2liouv(H,'right')`, `L_comm=hilb2liouv(H,'comm')`, `L_acomm=hilb2liouv(H,'acomm')`, and `R_vec=hilb2liouv(R,'statevec')`. The actions are compared respectively with `H*R`, `R*H`, `H*R-R*H`, and `H*R+R*H`, represented as column vectors. Each residual is a 2-norm; the maximum must be at most `1e-10`.
- For multiplicities `[2 3 5 8]`, compares `stevens(mult,1,0)`, `stevens(mult,1,1)`, and `stevens(mult,1,-1)` with `pauli(mult)`'s `z`, `x`, and `y` operators, respectively. The maximum of the three Frobenius-norm differences must be at most `1e-10`.

## Assumptions and limits

The Liouville identities are tested on one random matrix pair at dimension 6 per invocation, with no random seed set in the function. The rank-1 comparison covers only the four listed multiplicities and the stated component correspondence. The source specifies conditional success messages and tolerances but embeds no fixed numerical residuals; this entry does not assert that a run passed.

## Source

[`examples/fundamentals/operator_tests/commutation_8.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/operator_tests/commutation_8.m)
