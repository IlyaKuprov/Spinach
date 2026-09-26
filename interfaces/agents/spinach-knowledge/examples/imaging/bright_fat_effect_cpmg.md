# examples/imaging/bright_fat_effect_cpmg.m

- Signature: `bright_fat_effect_cpmg()`

## Purpose

Bright fat effect under CPMG echo train -magnetisation losses are greater in MRI experiments on J-coupled systems because co- herences are lost in the depths of the Hilbert space. Simulation time: minutes, faster with a Tesla V100 GPU.

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Bright fat effect under CPMG echo train -magnetisation losses
- are greater in MRI experiments on J-coupled systems because co-
- herences are lost in the depths of the Hilbert space.
- Simulation time: minutes, faster with a Tesla V100 GPU.
- Magnetic induction
- Spin system
- Chemical shifts
- J-coupling
- Spins 1,2,3 are molecule A
- and 4,5,6 are molecule B
- Kinetic rate matrix (Hz)
- Basis set
