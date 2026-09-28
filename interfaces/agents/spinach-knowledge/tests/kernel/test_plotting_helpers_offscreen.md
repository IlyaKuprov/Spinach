# tests/kernel/test_plotting_helpers_offscreen.m

- Signature: `result=test_plotting_helpers_offscreen()`

## Purpose

Tests offscreen execution of Spinach plotting helpers under invisible figures, without relying on image comparison.

## Physical / mathematical content

The test uses deterministic one-, two-, and three-dimensional spectra, MRI image data, a signed volume, and a compact mesh with heterogeneous Voronoi cells.

## Numerical / algorithmic content

Checks plotted data and frequency-axis sizes, graphics object counts and properties, figure dimensions, Voronoi plotting arrays, concentration-plot connectivity, and cap and side-wall areas.

## Outputs

- `result` — regression test result with explanatory messages.

## Implementation structure

- Announce the test target and initialize the regression result.
- Force invisible figures during the test, restoring figure visibility during cleanup.
- Build a minimal spin-system structure used by the plotting routines.
- Exercise house-style figure helpers and one-, two-, and three-dimensional spectral plotting.
- Exercise MRI, volume, COMSOL mesh, and concentration plotting utilities.