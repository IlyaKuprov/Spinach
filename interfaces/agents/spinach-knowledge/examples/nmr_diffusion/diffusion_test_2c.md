# examples/nmr_diffusion/diffusion_test_2c.m

- Signature: `diffusion_test_2c()`

## Purpose

A standard diffusion equation solver with no spin dynamics present. Anisotropic diffusion with a pe- riodic boundary condition. Calculation time: minutes.

## Physical / mathematical content

- Diffusion examples. The dominant mathematics is diffusion or advection-diffusion PDE propagation, sometimes with additional spin phase accumulation under gradients.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- A standard diffusion equation solver with no spin
- dynamics present. Anisotropic diffusion with a pe-
- riodic boundary condition.
- Calculation time: minutes.
- Load the phantom
- Ghost spin
- No spin interactions
- Basis set
- Spinach housekeeping
- Sample geometry
- 2D flow parameters
- 2D diffusion tensor field
