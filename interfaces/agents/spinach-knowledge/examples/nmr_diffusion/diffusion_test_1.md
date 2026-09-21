# examples/nmr_diffusion/diffusion_test_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_diffusion/diffusion_test_1.m`
- Signature: `diffusion_test_1()`
- Total lines: 61

## Purpose

A standard diffusion equation solver with no spin dynamics present. Calculation time: seconds

## Physical / mathematical content

- Diffusion examples. The dominant mathematics is diffusion or advection-diffusion PDE propagation, sometimes with additional spin phase accumulation under gradients.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- A standard diffusion equation solver with no spin
- dynamics present.
- Calculation time: seconds
- Ghost spin
- No spin interactions
- Basis set
- Spinach housekeeping
- Sample geometry
- Diffusion and flow parameters
- Diffusion and flow generator
- Initial condition
- Timing parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `v2fplanck()`, `inflate()`, `evolution()`, `kfigure()`, `traj()`, `kylabel()`, `kxlabel()`, `pause()`.
