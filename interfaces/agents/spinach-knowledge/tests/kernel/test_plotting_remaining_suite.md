# tests/kernel/test_plotting_remaining_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_plotting_remaining_suite.m`
- Signature: `result=test_plotting_remaining_suite()`
- Total lines: 457

## Purpose

Tests remaining Spinach plotting helper gaps under offscreen graphics. Syntax: result=test_plotting_remaining_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file also defines local helper function(s): `local_test_axes()`, `local_test_arrays()`, `local_test_graphics()`, `local_test_tensors()`, `local_test_interactive_guards()`, `local_test_write_movie()`, `local_expect_error()`, `local_plot_system()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `get()`, `set()`, `onCleanup()`, `local_cleanup()`, `local_plot_system()`, `local_test_axes()`, `local_test_arrays()`, `local_test_graphics()`, `local_test_tensors()`, `local_test_interactive_guards()`, `ft_axis()`, `test_close()`, `sweep2ticks()`, `fft_freq_axis()`, `ifft_time_axis()`.
