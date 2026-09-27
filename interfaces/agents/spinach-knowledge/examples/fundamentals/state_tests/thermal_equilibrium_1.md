# examples/fundamentals/state_tests/thermal_equilibrium_1.m

- Signature: `thermal_equilibrium_1()`

## Purpose

Check equilibrium magnetisations for `E8`, `1H`, `14N`, and `15N` against level-population results in Spinach's Zeeman Hilbert, Zeeman Liouville, and spherical-tensor Liouville formalisms.

## Method

The four-spin test system is set to 14.1 T and 4.2 K. Its isotopes are `E8`, `1H`, `14N`, and `15N`, with Zeeman values `{2.002319, 1.0, 2.0, 3.0}`. The script sets scalar couplings of `1e6` for pairs 1–2 and 2–3, `1e3` for pairs 1–3 and 3–4, and `1e2` for pair 1–4. Source comments say to remove these couplings to obtain a machine-precision match. For each formalism, the script obtains `equilibrium(spin_system)` and evaluates the four Lz expectation values. It compares these with analytical values computed using `levelpop`; any absolute discrepancy above `1e-5` fails the cross-formalism test.
