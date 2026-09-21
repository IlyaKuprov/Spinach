# tests/kernel/test_dynamic_examples_smoke.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_examples_smoke.m`
- Signature: `result=test_dynamic_examples_smoke()`
- Total lines: 149

## Purpose

Tests compact dynamic example-stage execution with plotting. Syntax: result=test_dynamic_examples_smoke()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The file also defines local helper function(s): `local_test_acquire_1d()`, `local_test_ct_cosy_2d()`, `local_cleanup()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `get()`, `set()`, `onCleanup()`, `local_cleanup()`, `local_test_acquire_1d()`, `local_test_ct_cosy_2d()`, `test_spin_system()`, `state()`, `liquid()`, `fftshift()`, `spectrum_ref()`, `fid_ref()`, `test_close()`, `kfigure()`, `plot_1d()`.
