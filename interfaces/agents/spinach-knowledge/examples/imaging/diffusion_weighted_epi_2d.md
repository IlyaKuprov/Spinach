# examples/imaging/diffusion_weighted_epi_2d.m

- Signature: `diffusion_weighted_epi_2d()`

## Purpose

2D echo planar imaging example in the presence of istropic diffusion. Stejskal-Tanner SE echo planar diffusion-weighted pulse sequence from Figure 1 in (https://doi.org/10.1148/radiol.09090021). Simulation time: minutes, faster with a Tesla V100 GPU.

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- 2D echo planar imaging example in the presence of istropic diffusion.
- Stejskal-Tanner SE echo planar diffusion-weighted pulse sequence from
- Figure 1 in (https://doi.org/10.1148/radiol.09090021).
- Simulation time: minutes, faster with a Tesla V100 GPU.
- Isotopes
- Magnetic induction
- Chemical shifts
- Relaxation model
- Disable path tracing
- This needs a GPU
- sys.enable={'gpu'};
- Basis set
