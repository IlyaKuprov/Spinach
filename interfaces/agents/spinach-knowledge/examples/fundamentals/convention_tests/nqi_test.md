# examples/fundamentals/convention_tests/nqi_test.m

- MATLAB implementation: [examples/fundamentals/convention_tests/nqi_test.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/convention_tests/nqi_test.m)

- Signature: `nqi_test()`

## Purpose

Exercises the reverse decomposition of a spin-1 Hamiltonian: `ham2nqi` extracts `omega` and `Q` from a randomly generated traceless Hermitian 3×3 matrix, then a Spinach ¹⁴N model is used to reconstruct that matrix. The coupling scaling and operator components are specified in the source; no additional unit interpretation is assumed here.

## Setup and reconstruction

Run `nqi_test()`. It creates `H_T` from a complex random 3×3 matrix, adds its conjugate transpose, and removes the trace with `remtrace`. It calls `[omega,Q]=ham2nqi(H_T)`, then builds a system with zero magnet field, isotope `14N`, and coupling matrix `Q/(2*pi)`. The basis uses `zeeman-hilb` formalism with `approximation='none'`; after `create` and `basis`, the code sets the lab frame, obtains `[H_iso,H_aniso]=hamiltonian(spin_system)`, and forms `H_S` from `H_iso + orientation(H_aniso,[0 0 0])` plus the coefficients `omega(1)`, `omega(2)`, and `omega(3)` multiplying the `Lx`, `Ly`, and `Lz` operators for `14N`.

The reconstruction passes when `norm(H_T-H_S,2)` is no greater than `1e-6*norm(H_T+H_S,2)`. This source-defined relative comparison is the acceptance criterion; it is not a reported fit outcome.

## Observable result and scope

The function prints `Quadrupolar reconstruction test PASSED.` when the comparison passes, or raises `Quadrupolar reconstruction test FAILED.` otherwise. It produces no plot and declares no output arguments. Each run uses a newly sampled random matrix; the source sets no seed and tests a single 3×3 case per invocation.
