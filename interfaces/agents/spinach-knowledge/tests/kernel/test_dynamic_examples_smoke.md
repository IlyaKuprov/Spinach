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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Dynamic plotting example smoke paths\n')`.
- Lines 20-23: State the dynamic example-stage target of the test; implemented by `result=new_test_result('examples/dynamic_examples_smoke', 'Dynamic plotting example smoke paths', 'compact example-stage calculations must run, process, and plot determi…`.
- Lines 25-26: Force invisible figures during the test; implemented by `old_visibility=get(groot,'defaultFigureVisible')`.
- Lines 30-31: Run and plot a compact one-dimensional acquisition path; implemented by `result=local_test_acquire_1d(result)`.
- Lines 33-34: Run and plot a compact two-dimensional CT-COSY path; implemented by `result=local_test_ct_cosy_2d(result)`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('examples/dynamic_examples_smoke', 'Dynamic plotting example smoke paths', 'compact example-stage calculations must run, process, and plot determi…`.
- Lines 26: computes `old_visibility` using `old_visibility=get(groot,'defaultFigureVisible')`.
- Lines 28: computes `cleaner` using `cleaner=onCleanup(@()local_cleanup(old_visibility))`.

### Local helper functions

- Line 39: `local_test_acquire_1d()` — `function result=local_test_acquire_1d(result)`. Build a zero-offset one-spin Liouville-space system
  - Representative operation: `sys.magnet=14.1`.
  - Representative operation: `sys.isotopes={'1H'}`.
- Line 86: `local_test_ct_cosy_2d()` — `function result=local_test_ct_cosy_2d(result)`. Build the two-spin system used in the CT-COSY plotting example
  - Representative operation: `sys.isotopes={'1H','1H'}`.
  - Representative operation: `sys.magnet=5.9`.
- Line 142: `local_cleanup()` — `function local_cleanup(old_visibility)`. Restore figure state after success or failure
  - Representative operation: `close all force`.
  - Representative operation: `set(groot,'defaultFigureVisible',old_visibility)`.

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
