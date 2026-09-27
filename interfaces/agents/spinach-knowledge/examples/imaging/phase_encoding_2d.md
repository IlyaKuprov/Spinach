# examples/imaging/phase_encoding_2d.m

- Signature: `phase_encoding_2d()`

## Purpose

Simple phase-encoded 2D imaging example. Calculation time: seconds. Ahmed Allami and Ilya Kuprov.

## Physical / mathematical content

- Simulates a single `1H` spin at 5.9 T with zero chemical shift and a `t1_t2` relaxation model. The relaxation rates are `R1 = 30.0` and `R2 = 70.0`, with diagonal relaxation terms retained and zero equilibrium.
- Loads `R1Ph` and `R2Ph` from `../../etc/phantoms/letter_a.mat` as spatial relaxation phantoms. The initial state is `Lz`, and the detection state is `L+`, each with a uniform spatial phantom.

## Numerical / algorithmic content

- Uses the `sphten-liouv` formalism with no basis approximation; path tracing and Krylov methods are disabled.
- Sets an image size of `[101 105]` and a spatial grid of `[108 90]` over dimensions `[0.30 0.25]`, with `{'period',3}` spatial derivatives. The readout gradient is `4.3e-3` T/m for `2e-3` s, the phase-encoding gradient is `3.8e-3` T/m for `1e-3` s, and the echo time is `0.025` s.

## Implementation structure

- Creates the spin system and basis, constructs relaxation operators with `rlx_t1_t2`, and runs `imaging(spin_system,@phase_enc_2d,parameters)`.
- Plots the recorded image alongside the `R1` and `R2` phantoms.
