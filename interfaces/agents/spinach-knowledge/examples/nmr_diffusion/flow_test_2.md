# examples/nmr_diffusion/flow_test_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_diffusion/flow_test_2.m`
- Signature: `flow_test_2()`
- Total lines: 62

## Purpose

A combination of diffusion and flow in two dimensions with a periodic boundary condition. Calculation time: minutes.

## Physical / mathematical content

- Diffusion examples. The dominant mathematics is diffusion or advection-diffusion PDE propagation, sometimes with additional spin phase accumulation under gradients.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- A combination of diffusion and flow in two dimensions
- with a periodic boundary condition.
- Calculation time: minutes.
- Load the phantom
- Ghost spin
- No spin interactions
- Basis set
- Spinach housekeeping
- Sample geometry
- 2D flow field
- 2D diffusion tensor field
- Diffusion and flow generator

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `create()`, `basis()`, `v2fplanck()`, `inflate()`, `evolution()`, `kfigure()`, `traj()`, `pause()`.
