# examples/fundamentals/state_tests/thermal_equilibrium_1.m

Source: [examples/fundamentals/state_tests/thermal_equilibrium_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/state_tests/thermal_equilibrium_1.m)

Signature: `thermal_equilibrium_1()`

## Tested question

For one four-spin configuration, do the four equilibrium Lz expectation values obtained with each of Spinach's three basis formalisms agree with single-spin level-population references within a stated absolute tolerance?

## System and reference

The source sets `sys.magnet=14.1`, `inter.temperature=4.2`, isotopes `E8`, `1H`, `14N`, and `15N`, and scalar Zeeman entries `{2.002319, 1.0, 2.0, 3.0}`. Its scalar couplings are 1–2 = `1e6`, 2–3 = `1e6`, 1–3 = `1e3`, 3–4 = `1e3`, and 1–4 = `1e2`. The source comments say to remove these couplings to obtain machine-precision agreement; it does not report a measured precision or numerical magnetisations.

For each of `zeeman-hilb`, `zeeman-liouv`, and `sphten-liouv`, with approximation `none`, the source forms `rho=equilibrium(spin_system)` and evaluates `trace(state(spin_system,'Lz',isotope)'*rho)` for each isotope. The reference vector comes from `levelpop`: the E8 populations are weighted by projections +3.5 through -3.5, the spin-1/2 populations for 1H and 15N by +0.5 and -0.5, and the 14N reference by `P(1)-P(3)`.

## Check, output, and limits

The four-by-three numerical array is compared with the four-by-one reference for every formalism. Any absolute difference greater than `1e-5` raises `Cross-formalism thermal equilibrium test FAILED.`; otherwise the function prints the corresponding PASSED message. It reports no individual values. This is an equilibrium and state-space comparison for the specified configuration, not a quadrature test or a general validation of interacting systems; the source's coupling-removal comment limits interpretation of the stated machine-precision match.
