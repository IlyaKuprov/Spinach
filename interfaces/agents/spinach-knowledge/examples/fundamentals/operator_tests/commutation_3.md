# examples/fundamentals/operator_tests/commutation_3.m

- Signature: `commutation_3()`
- Source: [examples/fundamentals/operator_tests/commutation_3.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/operator_tests/commutation_3.m)

## Purpose

Checks angular-momentum commutators for uranium and two-spin raising/lowering products in three Spinach formalisms.

## Model and operator identities

The no-argument function specifies a two-spin system with isotopes `1H` and `235U`, `sys.magnet=14.1`, scalar Zeeman entries `2.5` and `1.0`, and symmetric scalar coupling entries of `10`. It constructs a basis with `approximation='none'` for each of `zeeman-hilb`, `zeeman-liouv`, and `sphten-liouv`. For `235U`, the single-spin operators `L+`, `L-`, `Lx`, `Ly`, and `Lz` are tested against `[Lz,L+]=L+`, `[Lz,L-]=-L-`, and `[Lx,Ly]=i Lz`. The function also forms the two-spin products `L+_H L+_U` and `L-_H L-_U` in system-spin order; their commutators with the uranium and hydrogen `Lz` operators are checked with the corresponding positive and negative ladder signs.

## Calling and numerical checks

Call `commutation_3()` in MATLAB with the Spinach functions used by the source available. There are seven Frobenius-norm residuals per formalism, stored in a `7-by-3` array. The source's success branch prints `Cross-formalism commutation test PASSED.` when the Frobenius norm of the full array is below `1e-6`; otherwise it raises `Cross-formalism commutation test FAILED.` The script prints no individual residual values.

This is a fixed two-spin algebra check over the three listed representations, not a dynamics simulation or a general test over spin systems, coupling models, or approximations.
