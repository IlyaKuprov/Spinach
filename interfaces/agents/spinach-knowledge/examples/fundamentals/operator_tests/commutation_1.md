# examples/fundamentals/operator_tests/commutation_1.m

- MATLAB implementation: [examples/fundamentals/operator_tests/commutation_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/operator_tests/commutation_1.m)

- Signature: `commutation_1()`

## Purpose

Checks the angular-momentum operator commutation identities in three Spinach representations. This is an algebra and representation-convention check, not a magnetic-resonance propagation or experimental model.

## System and operators

The source creates one `1H` spin with `sys.magnet = 0` and zero scalar Zeeman interaction. For each formalism `zeeman-hilb`, `zeeman-liouv`, and `sphten-liouv`, it uses basis approximation `none` and obtains `Lx`, `Ly`, `Lz`, `L+`, and `L-` with `operator`.

The three residuals test `Lz*L+ - L+*Lz = L+`, `Lz*L- - L-*Lz = -L-`, and `Lx*Ly - Ly*Lx = 1i*Lz`. For each formalism, their Frobenius norms fill a column of a complex 3-by-3 answer array.

## Source-stated check and limits

The script prints `Cross-formalism commutation test PASSED.` only when `norm(answer,'fro') < 1e-6`; otherwise it raises an error reporting failure. This describes the source's conditional check, not an observed run: no MATLAB execution is claimed here. The test covers these three identities for this one-spin setup and these three formalisms; it does not test dynamics, other basis approximations, or broader operator behaviour.
