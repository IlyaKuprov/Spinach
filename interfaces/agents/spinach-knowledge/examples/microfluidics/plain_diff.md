# examples/microfluidics/plain_diff.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/microfluidics/plain_diff.m`
- Signature: `plain_diff()`
- Total lines: 115

## Purpose

Simple diffusion simulation without spin dynamics. Longitudinal magnetisation is tracked as a function of time.

## Physical / mathematical content

- Microfluidics examples. The coupled model is spin dynamics plus advection-diffusion-reaction transport on a mesh or regular grid. Numerical issues include finite-difference operators, mesh interpolation, and coupled reaction-flow evolution.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Simple diffusion simulation without spin dynamics. Longitudinal
- magnetisation is tracked as a function of time.
- Import hydrodynamics information
- One proton
- Chemical shift (water)
- Basis set
- Algorithmic switches
- Spinach housekeeping
- Initial condition: Lz in one cell in the middle
- Detection state: Lz in all cells
- Sequence and timing parameters
- Set assumptions

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `comsol_import()`, `create()`, `basis()`, `state()`, `assume()`, `hamiltonian()`, `frqoffset()`, `relaxation()`, `drainage()`, `speye()`, `meshflow()`, `fpl2phan()`, `traj()`, `kfigure()`, `scale_figure()`, `camproj()`.
