# examples/nmr_spen/ufcosy_3spin.m

- Signature: `ufcosy_3spin()`

## Purpose

Ultrafast COSY for a coupled three-spin system. Calculation time: minutes on NVidia Tesla A100, much longer on CPU Jean-Nicolas Dumez Ludmilla Guduff

## Physical / mathematical content

- The source defines a coupled three-spin system and sets up a 15 mm sample with 500 spatial points; diffusion and flow are set to zero.
- The imaging simulation uses the `@spencosy` sequence with encoding and coherence-selection parameters specified in the source.

## Numerical / algorithmic content

- The script Fourier-transforms the FID along its second dimension, then displays the magnitude as a contour plot.

## Implementation structure

- Constructs the spin system and basis, configures sample and sequence parameters, runs imaging, and plots the processed spectrum.
