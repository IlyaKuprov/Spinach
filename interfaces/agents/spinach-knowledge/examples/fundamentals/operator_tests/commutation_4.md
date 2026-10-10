# examples/fundamentals/operator_tests/commutation_4.m

- Signature: `commutation_4()`
- Source: [examples/fundamentals/operator_tests/commutation_4.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/operator_tests/commutation_4.m)

## Purpose

Checks the central-transition (CT) angular-momentum commutators in three Spinach formalisms.

## Model and operator identities

The no-argument function specifies `1H` and `235U`, with `sys.magnet=14.1`, scalar Zeeman entries `2.5` and `1.0`, and symmetric scalar coupling entries of `10`. It builds a basis with `approximation='none'` under `zeeman-hilb`, `zeeman-liouv`, and `sphten-liouv`. In each representation it obtains `CTx`, `CTy`, `CTz`, `CT+`, and `CT-` for `235U`, then checks `[CTz,CT+]=CT+`, `[CTz,CT-]=-CT-`, and `[CTx,CTy]=i CTz`. Although the two-spin system is built, the tested operators are the uranium CT operators; this example does not test mixed-spin products.

## Calling and numerical checks

Call `commutation_4()` in MATLAB with the Spinach functions used by the source available. The three Frobenius-norm residuals are stored for each formalism. The source's success branch prints `Cross-formalism commutation test PASSED.` when the norm of the full `3-by-3` residual array is below `1e-6`; otherwise it raises `Cross-formalism commutation test FAILED.`. No individual residuals are displayed.

The check covers these three CT identities, one uranium isotope, the specified two-spin setup, and the three listed representations; it is not a propagation test or a broad comparison of CT conventions.
