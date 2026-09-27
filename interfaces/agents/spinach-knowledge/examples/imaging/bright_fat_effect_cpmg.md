# examples/imaging/bright_fat_effect_cpmg.m

- Signature: `bright_fat_effect_cpmg()`

## Purpose

Simulates the bright fat effect under a CPMG echo train: magnetisation losses are greater in MRI experiments on J-coupled systems because coherences are lost in the depths of the Hilbert space. Simulation time is minutes, faster with a Tesla V100 GPU.

## Physical / mathematical content

At a magnetic induction of 3.0, the model has six 1H spins: spins 1–3 form molecule A and spins 4–6 form molecule B. Both have chemical shifts of 1.0, 2.0, and 3.0. The J-couplings within A are zero; those within B are 11, 17, and 23 for pairs 4–5, 4–6, and 5–6, respectively. The kinetic rate matrix is [0 0; 0 0] Hz, with concentrations [1 1]. Initial-state phantoms come from bright_fat_left.mat and bright_fat_right.mat; there is no diffusion or flow.

## Numerical / algorithmic content

Uses the sphten-liouv formalism with no basis approximation and disables path tracing. The CPMG simulation uses 48 pulses, a dec_time of 80e-3, zero offset, and no decoupling. The sample has dimensions [0.30 0.25], a [100 200] grid, and derivative setting {'period',3}. Relaxation phantoms and operators are empty. Imaging is run with imaging(spin_system,@cpmg_dec,parameters).

## Implementation structure

bright_fat_effect_cpmg() defines the spin system and basis, loads the two phantoms, and assigns initial Lz states to spins [1 2 3] and [4 5 6]. A uniform coil phantom detects the 1H Lx state. The function plots surf(abs(mri)) with a reversed X axis and the title 'Bright fat effect under CPMG echo train'. GPU enabling is present but commented out.
