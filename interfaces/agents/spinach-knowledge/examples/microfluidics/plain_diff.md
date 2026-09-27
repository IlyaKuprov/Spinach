# examples/microfluidics/plain_diff.m

- Signature: `plain_diff()`

## Purpose

Simple diffusion simulation without spin dynamics. Longitudinal magnetisation is tracked as a function of time.

## Physical / mathematical content

- This example isolates diffusion and distal-pipe drainage on an imported COMSOL mesh: it zeros both mesh velocity components, initializes longitudinal magnetisation in cells 1240 and 1246, and detects it with a uniform Lz coil. No spin-subspace dynamics or chemical reaction is included.

## Numerical / algorithmic content

- `meshflow` propagates the spatial signal with diffusion coefficient 1e-7 m^2/s and a drainage term in the distal cells; the script then plots the concentration-like trajectory over the mesh.

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
