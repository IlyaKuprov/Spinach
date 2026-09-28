# examples/imaging/echo_planar_2d.m

- Signature: `echo_planar_2d()`

## Purpose

Echo planar imaging example in 2D for a brain phantom. Simulation time: seconds, faster with a Tesla V100 GPU.

## Physical / mathematical content

- Simulates a `1H` spin system at a magnetic induction of 5.9 T, with a zero chemical shift and a `t1_t2` relaxation model. The initial longitudinal-magnetization and receive-coil transverse-magnetization states use `Lz` and `L+`, respectively.
- Uses one 2D slice (index 50) of the `brain-medres` R1, R2, and proton-density phantoms. The R1 and R2 maps supply spatially varying relaxation; the proton-density map supplies the initial-state phantom. Diffusion and flow are set to zero.

## Numerical / algorithmic content

- Creates the spin system in the `sphten-liouv` formalism without basis approximation, obtains the T1/T2 relaxation superoperators, and runs `imaging(spin_system,@epi_2d,parameters)`.
- Sets image size to `[101 105]`, readout and phase-encoding gradient amplitudes to `5.3e-3` and `4.8e-3` T/m, and both gradient durations to `2e-3` s. Spatial derivatives use `{'period',3}`.
- Halves the phase-encoding gradient amplitude after simulation for the field-of-view calculation, then plots the recorded image alongside the R1 and R2 phantom slices.

## Implementation structure

- Defines isotopes, field, chemical shift, relaxation, and basis settings; disables path tracing. The GPU enable line is present but commented out.
- Loads the phantom maps, configures geometry, relaxation and state phantoms, runs the EPI sequence, and plots its result.
