# examples/relaxation_theory/sat_rec_1.m

[Source file](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/sat_rec_1.m) · Signature: `sat_rec_1()`

## Purpose

A single-proton saturation-recovery example (calculation time: seconds). Its implemented relaxation model is phenomenological `t1_t2`, not a Bloch–Redfield or SLE calculation.

## Spin system and relaxation

The source sets `sys.magnet=14.1`, one `1H`, and Zeeman scalar `1.5`. It assigns both `r1_rates` and `r2_rates` the value `5.0`, uses the `dibari` equilibrium prescription, secular retention, and temperature value `298`. Units for these scalar inputs are not stated in the source. The basis is the complete `sphten-liouv` basis with no approximation.

## Preparation, evolution, and plotted signal

The initial state is `E` on `1H`; the detection coil is `Lz` on `1H`. The script constructs the static NMR Hamiltonian and relaxation superoperator, obtains an `Lx` pulse operator, then calls `step(...,pi)` on the initial state (the source comment calls this an inversion pulse). It evolves with `H+1i*R` for `1000` steps of `1e-3 s`, requesting the observable. The real answer is plotted against a `0–1 s` axis and labelled as the `S_Z` expectation value. The page reports the actual initial state and pulse operation without interpreting an uncomputed result.
