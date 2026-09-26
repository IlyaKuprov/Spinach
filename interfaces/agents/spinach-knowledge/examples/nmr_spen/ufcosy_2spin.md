# examples/nmr_spen/ufcosy_2spin.m

- Signature: `ufcosy_2spin()`

## Purpose

Ultrafast COSY for a coupled two-spin system. Calculation time: minutes on NVidia Tesla A100, much longer on CPU Jean-Nicolas Dumez Ludmilla Guduff

## Physical / mathematical content

- SPEN / ultrafast NMR examples. These files encode spatially dependent phase evolution and acquisition, linking pulse gradients, diffusion attenuation, and single-scan multidimensional encoding.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Implementation structure

- Ultrafast COSY for a coupled two-spin system.
- Calculation time: minutes on NVidia Tesla A100, much longer on CPU
- Jean-Nicolas Dumez
- Ludmilla Guduff
- Interactions
- Basis set
- Algorithmic options
- Spinach housekeeping
- Sample geometry
- Relaxation phantom
- Initial and detection state phantoms
- Diffusion and flow
