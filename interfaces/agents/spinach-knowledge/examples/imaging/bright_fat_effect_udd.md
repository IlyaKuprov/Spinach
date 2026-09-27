# examples/imaging/bright_fat_effect_udd.m

- Signature: `bright_fat_effect_udd()`

## Purpose

Simulates the bright fat effect under a UDD echo train: in MRI experiments on J-coupled systems, magnetisation losses are greater because coherences are lost in the depths of Hilbert space. Simulation time is minutes, faster with a Tesla V100 GPU.

## Physical / mathematical content

Models two three-spin 1H molecules at magnetic induction 3.0. Both have chemical shifts {1.0, 2.0, 3.0}; molecule A (spins 1, 2, 3) has zero pairwise J-couplings, while molecule B (spins 4, 5, 6) has couplings of 11, 17 and 23. The kinetic rate matrix is [0 0; 0 0] Hz, with concentrations [1 1]. Spatial initial-state phantoms are 1−left and 1−right, paired with the molecules’ Lz states; detection uses a uniform coil phantom and the 1H Lx state. Flow and diffusion are zero.

## Numerical / algorithmic content

Uses the sphten-liouv formalism without basis approximation. The UDD sequence has 48 pulses and a decoupling time of 80e-3, with zero offset and no decoupled spins. The spatial grid has dimensions [0.30 0.25], [100 200] points and derivative settings {'period',3}. No relaxation phantoms or operators are supplied. The simulation calls imaging(spin_system,@udd_dec,parameters).

## Implementation structure

Creates the spin system and basis, disables path tracing, and leaves GPU enablement commented out. Loads bright_fat_left.mat and bright_fat_right.mat, configures the sequence and spatial phantoms, then plots abs(mri) as a surface with the X direction reversed and the title 'Bright fat effect under UDD echo train'.
