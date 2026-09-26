# tests/kernel/test_finite_difference_suite.m

- Signature: `result=test_finite_difference_suite()`

## Purpose

Tests finite-difference and spectral differentiation helpers. Syntax: result=test_finite_difference_suite()

## Physical / mathematical content

- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.

## Outputs

- result -regression test result with explanatory messages
- The test checks finite-difference weights, finite-difference matrices,
- Fourier differentiation, Laplacians, FFT differentiation kernels,
- pseudomodulation, and matrix-exponential directional derivatives
- against exact simple cases.

## Implementation structure

- Tests finite-difference and spectral differentiation helpers. Syntax:
- result=test_finite_difference_suite()
- result -regression test result with explanatory messages
- The test checks finite-difference weights, finite-difference matrices,
- Fourier differentiation, Laplacians, FFT differentiation kernels,
- pseudomodulation, and matrix-exponential directional derivatives
- against exact simple cases.
- Announce the test target
- State the differentiation target of the test
- Three-point centred finite-difference weights at zero are exact and familiar
- Five-point wall finite-difference matrix differentiates quadratics exactly on a unit grid
- Savitzky-Golay differentiation recovers a cubic exactly on a uniform grid
