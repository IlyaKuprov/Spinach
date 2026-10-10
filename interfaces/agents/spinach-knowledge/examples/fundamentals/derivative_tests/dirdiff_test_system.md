# examples/fundamentals/derivative_tests/dirdiff_test_system.m

- MATLAB implementation: [examples/fundamentals/derivative_tests/dirdiff_test_system.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/derivative_tests/dirdiff_test_system.m)

- Signature: `[spin_system,Sx,Sy,Sz,Lx,Ly,H]=dirdiff_test_system(formalism)`

## Purpose and callable context

This helper constructs the spin system, normalised states, control operators, and drift Hamiltonian used by the directional-derivative examples. Call it with one of the three accepted character-vector values `sphten-liouv`, `zeeman-liouv`, or `zeeman-hilb`; any other input fails the source's consistency check. `dirdiff_8_rect()` is one mapped caller.

## System and basis

The source uses 100 non-interacting 13C spins for `sphten-liouv`, and 2 spins for either Zeeman formalism. The field is set to `28.18` T. All spins are 13C, with scalar chemical shifts equally spaced from `-100` to `+100` ppm; no explicit spin-spin couplings are assigned in this helper. Output is suppressed with `sys.output='hush'`.

For `sphten-liouv`, the basis uses approximation `IK-2`, proximity level 1, and connectivity `scalar_couplings`. For `zeeman-liouv` and `zeeman-hilb`, the basis uses approximation `none`. The function constructs the system with `create` and applies the selected basis with `basis`.

## Returned operators and states

`Sx`, `Sy`, and `Sz` are obtained from `state` using `Lx`, `Ly`, and `Lz` for 13C, then each is divided by `norm(full(state),2)`. The returned `Lx` and `Ly` are 13C operators from `operator`; `H` is `hamiltonian(assume(spin_system,'nmr'))`. These are generated model objects, not numerical derivative-test results. The helper contains no numerical comparison or pass/fail output; its explicit validation is limited to the formalism input check. Source: [examples/fundamentals/derivative_tests/dirdiff_test_system.m](../../../../../../examples/fundamentals/derivative_tests/dirdiff_test_system.m).