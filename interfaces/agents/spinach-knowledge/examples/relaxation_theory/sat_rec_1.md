# examples/relaxation_theory/sat_rec_1.m

- Signature: `sat_rec_1()`

## Purpose

Demonstrates a saturation-recovery experiment for a single proton, using a pi rotation to invert the initial state and then tracking longitudinal magnetization during relaxation. Calculation time: seconds.

## Physical / mathematical content

The single-`1H` spin is at `14.1 T`, with Zeeman scalar `1.5`. The model uses phenomenological `T1/T2` relaxation rates of `5.0 s^-1` each, the Di Bari equilibrium state, and temperature `298 K`.

## Numerical / algorithmic content

The script uses the complete `sphten-liouv` basis with no approximation. It creates the unit initial state and an `Lz` detection state, builds the static NMR Hamiltonian and relaxation superoperator, and applies a `pi` pulse about `Lx`. It then evolves under `H + 1i*R` for `1000` steps of `1e-3 s`, plots the real signal against a `0–1 s` time axis, and labels it as the `S_Z` expectation value.
