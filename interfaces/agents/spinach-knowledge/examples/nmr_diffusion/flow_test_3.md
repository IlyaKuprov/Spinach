# examples/nmr_diffusion/flow_test_3.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_diffusion/flow_test_3.m`
- Signature: `flow_test_3()`
- Total lines: 81

## Purpose

Circular flow in three-dimensional space in the absence of spin dynamics. Calculation time: minutes, faster on GPU.

## Physical / mathematical content

- Diffusion examples. The dominant mathematics is diffusion or advection-diffusion PDE propagation, sometimes with additional spin phase accumulation under gradients.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- Circular flow in three-dimensional space in the absence
- of spin dynamics.
- Calculation time: minutes, faster on GPU.
- Ghost spin
- No spin interactions
- Basis set
- Algorithmic options
- Spinach housekeeping
- Sample geometry
- Get a 3D grid
- Get circular wind vectors
- Constant diffusion tensor field

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `v2fplanck()`, `inflate()`, `evolution()`, `kfigure()`, `volplot()`, `traj()`.
