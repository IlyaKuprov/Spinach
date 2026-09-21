# examples/nmr_spen/ufdosy_1spin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_spen/ufdosy_1spin.m`
- Signature: `ufdosy_1spin()`
- Total lines: 87

## Purpose

Ultrafast DOSY for one spin. Calculation time: seconds on NVidia Tesla A100, much longer on CPU Ludmilla Guduff Jean-Nicolas Dumez

## Physical / mathematical content

- SPEN / ultrafast NMR examples. These files encode spatially dependent phase evolution and acquisition, linking pulse gradients, diffusion attenuation, and single-scan multidimensional encoding.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.

## Implementation structure

- Ultrafast DOSY for one spin.
- Calculation time: seconds on NVidia Tesla A100, much longer on CPU
- Ludmilla Guduff
- Jean-Nicolas Dumez
- Spin system
- Interactions
- Basis set
- Algorithmic options
- Spinach housekeeping
- Assumptions
- Sample geometry
- Relaxation phantom

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `assume()`, `state()`, `imaging()`, `fftshift()`, `spin()`, `kfigure()`, `kxlabel()`, `kylabel()`.
