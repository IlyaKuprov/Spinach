# tests/kernel/test_plotting_helpers_offscreen.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_plotting_helpers_offscreen.m`
- Signature: `result=test_plotting_helpers_offscreen()`
- Total lines: 404

## Purpose

Tests offscreen execution of Spinach plotting helpers. Syntax: result=test_plotting_helpers_offscreen()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file also defines local helper function(s): `local_test_house_style()`, `local_test_plot_1d()`, `local_test_plot_2d()`, `local_test_plot_3d()`, `local_test_misc_plots()`, `local_test_comsol_plots()`, `local_plot_system()`, `local_spectrum_2d()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `get()`, `set()`, `onCleanup()`, `local_cleanup()`, `local_plot_system()`, `local_test_house_style()`, `local_test_plot_1d()`, `local_test_plot_2d()`, `local_test_plot_3d()`, `local_test_misc_plots()`, `local_test_comsol_plots()`, `kfigure()`, `test_true()`, `ishandle()`, `strcmp()`.
