# examples/nmr_spen/ufdosy_1spin.m

- Signature: `ufdosy_1spin()`

## Purpose

Ultrafast DOSY for one spin. Calculation time: seconds on NVidia Tesla A100, much longer on CPU Ludmilla Guduff Jean-Nicolas Dumez

## Physical / mathematical content

- The source defines one 1H spin at 14.1 T with a 7.0 ppm shift. The sample is 15 mm long with 3000 spatial points; diffusion is set to 8e-10 m^2/s and flow to zero.
- The imaging simulation calls the `@spendosy` sequence and specifies the acquisition and encoding gradients in the source.

## Numerical / algorithmic content

- The script Fourier-transforms the FID along the spatial axis, then plots the magnitude against chemical shift and displacement coordinates derived from the acquisition and gradient parameters.

## Implementation structure

- Builds the spin system and basis, configures sample, diffusion, acquisition, and encoding parameters, runs imaging, and transforms and plots the result.
