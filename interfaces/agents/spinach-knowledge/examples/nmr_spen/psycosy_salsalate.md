# examples/nmr_spen/psycosy_salsalate.m

- Signature: `psycosy_salsalate()`

## Purpose

PSYCOSY of one salsalate ring. Calculation time: minutes on NVidia Tesla A100, much longer on CPU

## Physical / mathematical content

- SPEN / ultrafast NMR examples. These files encode spatially dependent phase evolution and acquisition, linking pulse gradients, diffusion attenuation, and single-scan multidimensional encoding.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- PSYCOSY of one salsalate ring.
- Calculation time: minutes on NVidia Tesla A100, much longer on CPU
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
