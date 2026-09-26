# tests/kernel/test_dynamic_examples_smoke.m

- Signature: `result=test_dynamic_examples_smoke()`

## Purpose

Tests compact dynamic example-stage execution with plotting. Syntax: result=test_dynamic_examples_smoke()

## Physical / mathematical content

- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Outputs

- result -regression test result with explanatory messages
- The test runs short liquid-state NMR calculations adapted from plotting
- examples, processes the deterministic signals, and verifies the plotted
- graphics objects under invisible offscreen figures.

## Implementation structure

- Tests compact dynamic example-stage execution with plotting. Syntax:
- result=test_dynamic_examples_smoke()
- result -regression test result with explanatory messages
- The test runs short liquid-state NMR calculations adapted from plotting
- examples, processes the deterministic signals, and verifies the plotted
- graphics objects under invisible offscreen figures.
- Announce the test target
- State the dynamic example-stage target of the test
- Force invisible figures during the test
- Run and plot a compact one-dimensional acquisition path
- Run and plot a compact two-dimensional CT-COSY path
- Build a zero-offset one-spin Liouville-space system
