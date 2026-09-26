# examples/imaging/fast_echo_2d.m

- Signature: `fast_echo_2d()`

## Purpose

Fast (in the experiment duration sense) spin echo 2D brain imaging example. Simulation time: hours, faster with a Tesla V100 GPU.

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Fast (in the experiment duration sense) spin echo 2D brain
- imaging example.
- Simulation time: hours, faster with a Tesla V100 GPU.
- Isotopes
- Magnetic induction
- Chemical shifts
- Relaxation model
- Disable path tracing
- This needs a GPU
- sys.enable={'gpu'};
- Basis set
- Spinach housekeeping
