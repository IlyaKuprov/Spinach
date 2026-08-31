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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Offscreen plotting helpers\n')`.
- Lines 20-23: State the plotting-helper target of the test; implemented by `result=new_test_result('kernel/plotting_helpers_offscreen', 'Offscreen plotting helpers', 'plotting helpers must create deterministic graphics objects under invisible of…`.
- Lines 25-26: Force invisible figures during the test; implemented by `old_visibility=get(groot,'defaultFigureVisible')`.
- Lines 30-31: Build a minimal spin-system structure used only by plotting routines; implemented by `spin_system=local_plot_system()`.
- Lines 33-34: Exercise house-style figure helpers; implemented by `result=local_test_house_style(result)`.
- Lines 36-37: Exercise one-dimensional spectral plotting; implemented by `result=local_test_plot_1d(result,spin_system)`.
- Lines 39-40: Exercise two-dimensional contour and stack plotting; implemented by `result=local_test_plot_2d(result,spin_system)`.
- Lines 42-43: Exercise three-dimensional spectrum plotting; implemented by `result=local_test_plot_3d(result,spin_system)`.
- Lines 45-46: Exercise MRI and volume plotting utilities; implemented by `result=local_test_misc_plots(result)`.
- Lines 48-49: Exercise COMSOL mesh and concentration plotting; implemented by `result=local_test_comsol_plots(result)`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/plotting_helpers_offscreen', 'Offscreen plotting helpers', 'plotting helpers must create deterministic graphics objects under invisible of…`.
- Lines 26: computes `old_visibility` using `old_visibility=get(groot,'defaultFigureVisible')`.
- Lines 28: computes `cleaner` using `cleaner=onCleanup(@()local_cleanup(old_visibility))`.
- Lines 31: computes `spin_system` using `spin_system=local_plot_system()`.

### Local helper functions

- Line 54: `local_test_house_style()` — `function result=local_test_house_style(result)`. Create an invisible figure and check default handling
  - Representative operation: `fig=kfigure('Visible','off')`.
  - Representative operation: `result=test_true(result,'kfigure invisible handle',ishandle(fig)&&strcmp(get(fig,'Visible'),'off'), 'kfigure returns a valid invisible figure handle when requested')`.
- Line 97: `local_test_plot_1d()` — `function result=local_test_plot_1d(result,spin_system)`. Define a deterministic one-dimensional spectrum
  - Representative operation: `parameters.spins={'1H'}`.
  - Representative operation: `parameters.offset=0`.
- Line 139: `local_test_plot_2d()` — `function result=local_test_plot_2d(result,spin_system)`. Define a deterministic two-dimensional spectrum
  - Representative operation: `parameters.spins={'1H','13C'}`.
  - Representative operation: `parameters.offset=[0 0]`.
- Line 187: `local_test_plot_3d()` — `function result=local_test_plot_3d(result,spin_system)`. Define a compact deterministic three-dimensional spectrum
  - Representative operation: `parameters.spins={'1H','13C','15N'}`.
  - Representative operation: `parameters.offset=[0 0 0]`.
- Line 216: `local_test_misc_plots()` — `function result=local_test_misc_plots(result)`. Exercise MRI image and k-space plotting branches
  - Representative operation: `fig=kfigure('Visible','off')`.
  - Representative operation: `parameters.spins={'1H'}`.
- Line 251: `local_test_comsol_plots()` — `function result=local_test_comsol_plots(result)`. Build a compact mesh with heterogeneous Voronoi cells
  - Representative operation: `vertices=[0 0;1 0;0 1;1 0;2 0;2 1;1 1;3 0;3 1;2 1]`.
  - Representative operation: `cells={[1 2 3],[4 5 6 7],[8 9 10]}`.
- Line 376: `local_plot_system()` — `function spin_system=local_plot_system()`. Provide the minimal fields used by plotting utilities
  - Representative operation: `spin_system.comp.isotopes={'1H','13C','15N'}`.
  - Representative operation: `spin_system.inter.magnet=14.1`.
- Line 387: `local_spectrum_2d()` — `function spectrum=local_spectrum_2d(nrows,ncols)`. Build smooth positive and negative peaks without random numbers
  - Representative operation: `[row_grid,col_grid]=ndgrid(linspace(-1,1,nrows),linspace(-1,1,ncols))`.
  - Representative operation: `pos_peak=exp(-8*((row_grid-0.30).^2+(col_grid+0.20).^2))`.

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
