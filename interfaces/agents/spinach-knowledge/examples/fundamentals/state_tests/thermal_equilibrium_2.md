# examples/fundamentals/state_tests/thermal_equilibrium_2.m

Source: [examples/fundamentals/state_tests/thermal_equilibrium_2.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/state_tests/thermal_equilibrium_2.m)

Signature: `thermal_equilibrium_2()`

## Tested question

Does the finite-temperature equilibrium state for one four-spin system represent the same Hilbert-space matrix when calculated in Zeeman Hilbert, Zeeman Liouville, and spherical-tensor Liouville formalisms?

## System and method

The source sets a 14.1 T field, 4.2 K temperature, isotopes `E8`, `1H`, `14N`, and `15N`, and scalar Zeeman entries `{2.002319, 1.0, 2.0, 3.0}`. Scalar couplings are 1–2 = `1e6`, 2–3 = `1e6`, 1–3 = `1e3`, 3–4 = `1e3`, and 1–4 = `1e2`. The source comments that removing these couplings gives machine-precision agreement; no numerical residuals are printed.

With approximation `none`, it calculates `equilibrium(spin_system)` in each formalism. The Zeeman-Liouville state is reshaped to 96-by-96. The spherical-tensor state is transformed by `sphten2zeeman(spin_system)`, reshaped to 96-by-96, and divided by the product of the spin multiplicities. The direct Zeeman-Hilbert matrix is the reference representation.

## Check, output, and limits

The source evaluates the matrix 2-norms of the Hilbert-vs-Zeeman-Liouville and Zeeman-Liouville-vs-projected-spherical-tensor differences. Either norm above `1e-8` raises `Cross-formalism thermal equilibrium test FAILED.`; otherwise the function prints the corresponding PASSED message. This is one finite-temperature state-space comparison, not a quadrature test, and it does not provide independent level-population values or establish agreement beyond this specified system and the source's stated coupling caveat.
