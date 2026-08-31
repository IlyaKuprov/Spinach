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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Remaining offscreen plotting helpers\n')`.
- Lines 20-23: State the remaining plotting-helper target of the test; implemented by `result=new_test_result('kernel/plotting_remaining_suite', 'Remaining plotting helper gaps', 'remaining plotting helpers must return deterministic arrays, graphics object…`.
- Lines 25-26: Force invisible figures during the test; implemented by `old_visibility=get(groot,'defaultFigureVisible')`.
- Lines 31-32: Build a minimal spin-system structure used by plotting routines; implemented by `spin_system=local_plot_system()`.
- Lines 34-35: Exercise deterministic one-dimensional and transform axes; implemented by `result=local_test_axes(result,spin_system)`.
- Lines 37-38: Exercise colour maps, contour levels, cropping, and volume zooming; implemented by `result=local_test_arrays(result,spin_system)`.
- Lines 40-41: Exercise ordinary graphics helpers under invisible figures; implemented by `result=local_test_graphics(result,spin_system)`.
- Lines 43-44: Exercise tensor-display graphics under invisible figures; implemented by `result=local_test_tensors(result)`.
- Lines 46-47: Exercise non-interactive integration and guarded interactive paths; implemented by `result=local_test_interactive_guards(result,spin_system)`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/plotting_remaining_suite', 'Remaining plotting helper gaps', 'remaining plotting helpers must return deterministic arrays, graphics object…`.
- Lines 26: computes `old_visibility` using `old_visibility=get(groot,'defaultFigureVisible')`.
- Lines 27: computes `warning_state` using `warning_state=warning('query','MATLAB:griddedInterpolant:MeshgridEval2DWarnId')`.
- Lines 29: computes `cleaner` using `cleaner=onCleanup(@()local_cleanup(old_visibility,warning_state))`.
- Lines 32: computes `spin_system` using `spin_system=local_plot_system()`.

### Local helper functions

- Line 52: `local_test_axes()` — `function result=local_test_axes(result,spin_system)`. Check even-point Fourier axis folding
  - Representative operation: `axis_even=ft_axis(1,8,4)`.
  - Representative operation: `result=test_close(result,'ft_axis even points',axis_even,[-3 -1 1 3],1e-12,1e-12, 'even Fourier axes drop the folded edge frequency')`.
- Line 121: `local_test_arrays()` — `function result=local_test_arrays(result,spin_system)`. Check blue-white-red colour map anchors
  - Representative operation: `colour_map=bwr_cmap()`.
  - Representative operation: `result=test_true(result,'bwr_cmap size',isequal(size(colour_map),[255 3]), 'bwr_cmap returns a 255-by-3 MATLAB colour map')`.
- Line 206: `local_test_graphics()` — `function result=local_test_graphics(result,spin_system)`. Exercise molecular stick plotting with supplied connectivity
  - Representative operation: `fig=figure('Visible','off')`.
  - Representative operation: `xyz=[0 0 0;1 0 0;1 1 0]`.
- Line 254: `local_test_tensors()` — `function result=local_test_tensors(result)`. Build minimal tensor-display input structures
  - Representative operation: `props=local_tensor_props()`.
  - Representative operation: `options.style='ellipsoids'`.
- Line 300: `local_test_interactive_guards()` — `function result=local_test_interactive_guards(result,spin_system)`. Build deterministic data for non-interactive 2D integration
  - Representative operation: `parameters.spins={'1H','13C'}`.
  - Representative operation: `parameters.offset=[0 0]`.
- Line 348: `local_test_write_movie()` — `function result=local_test_write_movie(result)`. Prepare a tiny three-dimensional figure for movie capture
  - Representative operation: `fig=figure('Visible','off')`.
  - Representative operation: `plot3([0 1],[0 1],[0 1],'k-')`.
- Line 367: `local_expect_error()` — `function result=local_expect_error(result,label,call_handle,message_part,why)`. Run the supplied call and require a matching failure
  - Representative operation: `try`.
  - Representative operation: `call_handle()`.
- Line 385: `local_plot_system()` — `function spin_system=local_plot_system()`. Provide the minimal fields used by plotting utilities
  - Representative operation: `spin_system.comp.isotopes={'1H','13C'}`.
  - Representative operation: `spin_system.inter.magnet=14.1`.

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
