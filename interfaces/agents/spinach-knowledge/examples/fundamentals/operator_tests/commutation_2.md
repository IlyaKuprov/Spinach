# examples/fundamentals/operator_tests/commutation_2.m

- MATLAB implementation: [examples/fundamentals/operator_tests/commutation_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/operator_tests/commutation_2.m)

- Signature: `commutation_2()`

## Purpose

Checks the angular-momentum operator commutation identities for a `235U` spin in three Spinach representations. It is an algebra and representation-convention check, not a magnetic-resonance propagation or experimental model.

## System and operators

The source creates one `235U` isotope with `sys.magnet = 0` and zero scalar Zeeman interaction. For each formalism `zeeman-hilb`, `zeeman-liouv`, and `sphten-liouv`, it uses basis approximation `none` and obtains `L+`, `L-`, `Lx`, `Ly`, and `Lz` with `operator`.

The three residuals test `Lz*L+ - L+*Lz = L+`, `Lz*L- - L-*Lz = -L-`, and `Lx*Ly - Ly*Lx = 1i*Lz`. Their Frobenius norms fill a complex 3-by-3 answer array, with one column for each formalism.

## Source-stated check and limits

The script prints `Cross-formalism commutation test PASSED.` only when `norm(answer,'fro') < 1e-6`; otherwise it raises an error reporting failure. This is the source's conditional pass/fail rule, not an observed run: no MATLAB execution is claimed here. The scope is limited to these identities, this one-isotope system, and these three formalisms; no dynamics or broader operator validation is reported.
