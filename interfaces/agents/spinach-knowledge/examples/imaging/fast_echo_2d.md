# examples/imaging/fast_echo_2d.m

- Signature: `fast_echo_2d()`

## Purpose

Fast (in the experiment duration sense) spin echo 2D brain imaging example. Simulation time: hours, faster with a Tesla V100 GPU.

## Physical / mathematical content

- Simulates a `1H` spin system at a magnetic induction of 5.9 with zero chemical shift. The `t1_t2` relaxation model uses diagonal relaxation, zero equilibrium, and R1 and R2 rates of 1.
- Uses slice 50 of the `brain-medres` R1, R2, and proton-density phantoms. The R1 and R2 slices supply spatial relaxation maps; the proton-density slice supplies the initial-state map.

## Numerical / algorithmic content

- Builds a `sphten-liouv` basis with no approximation and disables path tracing. GPU enablement is suggested in a comment but is not active in the source.
- Runs `imaging(spin_system,@fse,parameters)` for a `101 × 105` image, with a 5.3e-3 T/m readout gradient lasting 2e-3 s and a 4.8e-3 T/m phase-encoding gradient lasting 1e-3 s. The sequence parameters specify `1H`, no decoupling, and zero offset.
- Uses the phantom's first two dimensions and point counts for the 2D sample geometry, with periodic third-order differentiation (`{'period',3}`). The initial spin state is `Lz`; the detection state is `L+` with a uniform coil map.

## Implementation structure

- Creates the spin system and basis, obtains R1 and R2 relaxation superoperators, and loads the brain phantom maps.
- Simulates the `fse` sequence and plots the recorded image alongside the R1 and R2 phantom slices.
