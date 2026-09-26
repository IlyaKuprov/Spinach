# tests/kernel/test_plotting_remaining_suite.m

- Signature: `result=test_plotting_remaining_suite()`

## Purpose

Tests remaining Spinach plotting helper gaps under offscreen graphics. Syntax: result=test_plotting_remaining_suite()

## Physical / mathematical content

- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test exercises deterministic axis, contour, cropping, molecular,
- tensor-display, ultrafast, and guarded interactive plotting helpers under
- invisible figures without relying on image comparison.

## Implementation structure

- Tests remaining Spinach plotting helper gaps under offscreen graphics. Syntax:
- result=test_plotting_remaining_suite()
- result -regression test result with explanatory messages
- The test exercises deterministic axis, contour, cropping, molecular,
- tensor-display, ultrafast, and guarded interactive plotting helpers under
- invisible figures without relying on image comparison.
- Announce the test target
- State the remaining plotting-helper target of the test
- Force invisible figures during the test
- Build a minimal spin-system structure used by plotting routines
- Exercise deterministic one-dimensional and transform axes
- Exercise colour maps, contour levels, cropping, and volume zooming
