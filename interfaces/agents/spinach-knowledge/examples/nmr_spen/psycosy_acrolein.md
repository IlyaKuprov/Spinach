# examples/nmr_spen/psycosy_acrolein.m

- Signature: `psycosy_acrolein()`

## Purpose

PSYCOSY of Acrolein. Calculation time: hours, faster on a GPU.

## Physical / mathematical content

- SPEN / ultrafast NMR examples. These files encode spatially dependent phase evolution and acquisition, linking pulse gradients, diffusion attenuation, and single-scan multidimensional encoding.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- PSYCOSY of Acrolein.
- Calculation time: hours, faster on a GPU.
- Magnet
- Spin system
- Interactions
- Algorithmic options
- Basis set
- Spinach housekeeping
- Sample geometry
- Diffusion and flow
- Relaxation phantom
- Initial and detection state phantoms
