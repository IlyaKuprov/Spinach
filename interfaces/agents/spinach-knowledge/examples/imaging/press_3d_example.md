# examples/imaging/press_3d_example.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/imaging/press_3d_example.m`
- Signature: `press_3d_example()`
- Total lines: 79

## Purpose

PRESS excitation profile in three dimensions with tilted gradient system. Change the frequency under Pulse Parame- ters to move the hot spot through the sample. Simulation time: hours, faster with a Tesla V100 GPU.

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- PRESS excitation profile in three dimensions with tilted
- gradient system. Change the frequency under Pulse Parame-
- ters to move the hot spot through the sample.
- Simulation time: hours, faster with a Tesla V100 GPU.
- Magnetic induction
- Spin systems
- Basis set
- Disable path tracing
- This is here essential
- Spinach housekeeping
- Sequence parameters
- Pulse parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `relaxation()`, `state()`, `imaging()`, `volplot()`, `kxlabel()`, `kylabel()`, `kzlabel()`.
