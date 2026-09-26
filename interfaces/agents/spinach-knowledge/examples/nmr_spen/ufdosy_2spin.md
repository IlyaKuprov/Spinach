# examples/nmr_spen/ufdosy_2spin.m

- Signature: `ufdosy_2spin()`

## Purpose

Ultrafast DOSY for two coupled spins with additional complications like DD and CSA relaxation, and spatial flow. Calculation time: minutes on NVidia Tesla A100, much longer on CPU Ludmilla Guduff Jean-Nicolas Dumez Ilya Kuprov

## Physical / mathematical content

- SPEN / ultrafast NMR examples. These files encode spatially dependent phase evolution and acquisition, linking pulse gradients, diffusion attenuation, and single-scan multidimensional encoding.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.

## Implementation structure

- Ultrafast DOSY for two coupled spins with additional complications
- like DD and CSA relaxation, and spatial flow.
- Calculation time: minutes on NVidia Tesla A100, much longer on CPU
- Ludmilla Guduff
- Jean-Nicolas Dumez
- Ilya Kuprov
- Spin system
- Interactions
- Basis set
- Relaxation theory
- Algorithmic options
- Spinach housekeeping
