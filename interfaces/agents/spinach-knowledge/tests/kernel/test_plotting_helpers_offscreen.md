# tests/kernel/test_plotting_helpers_offscreen.m

- Signature: `result=test_plotting_helpers_offscreen()`

## Purpose

Tests offscreen execution of Spinach plotting helpers. Syntax: result=test_plotting_helpers_offscreen()

## Physical / mathematical content

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Outputs

- result -regression test result with explanatory messages
- The test exercises plotting helpers under invisible figures, and checks
- graphics object creation, axis sizes, returned data arrays, and figure
- helper side effects without relying on image comparison.

## Implementation structure

- Tests offscreen execution of Spinach plotting helpers. Syntax:
- result=test_plotting_helpers_offscreen()
- result -regression test result with explanatory messages
- The test exercises plotting helpers under invisible figures, and checks
- graphics object creation, axis sizes, returned data arrays, and figure
- helper side effects without relying on image comparison.
- Announce the test target
- State the plotting-helper target of the test
- Force invisible figures during the test
- Build a minimal spin-system structure used only by plotting routines
- Exercise house-style figure helpers
- Exercise one-dimensional spectral plotting
