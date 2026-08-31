# tests/kernel/test_dynamic_grid_plot.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_grid_plot.m`
- Signature: `result=test_dynamic_grid_plot()`
- Total lines: 85

## Purpose

Tests grid_plot() under offscreen graphics. Syntax: result=test_dynamic_grid_plot()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file also defines local helper function(s): `local_patch_colours()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Spherical grid plotting\n')`.
- Lines 19-22: State the grid_plot target of the test; implemented by `result=new_test_result('kernel/dynamic_grid_plot', 'Offscreen spherical grid plotting', 'grid_plot() must draw one patch per Voronoi cell and honour centre-dot options.')`.
- Lines 24-25: Force invisible figures during the test; implemented by `old_visibility=get(groot,'defaultFigureVisible')`.
- Lines 29-30: Build a regular tetrahedral grid on the unit sphere; implemented by `xyz=[1 1 -1 -1;1 -1 1 -1;1 -1 -1 1]`.
- Lines 36-37: Compute the spherical Voronoi tessellation once; implemented by `[~,~,vorn]=voronoisphere(xyz)`.
- Lines 39-40: Draw supplied tessera with numeric colours and no centre dots; implemented by `fig=figure('Visible','off')`.
- Lines 55-56: Draw with internally generated tessera and default centre dots; implemented by `fig=figure('Visible','off')`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/dynamic_grid_plot', 'Offscreen spherical grid plotting', 'grid_plot() must draw one patch per Voronoi cell and honour centre-dot options.')`.
- Lines 25: computes `old_visibility` using `old_visibility=get(groot,'defaultFigureVisible')`.
- Lines 27: computes `cleaner` using `cleaner=onCleanup(@()set(groot,'defaultFigureVisible',old_visibility))`.
- Lines 30: computes `xyz` using `xyz=[1 1 -1 -1;1 -1 1 -1;1 -1 -1 1]`.
- Lines 32: computes `x` using `x=xyz(1,:).'`.
- Lines 33: computes `y` using `y=xyz(2,:).'`.
- Lines 34: computes `z` using `z=xyz(3,:).'`.
- Lines 37: computes `[~,~,vorn]` using `[~,~,vorn]=voronoisphere(xyz)`.
- Lines 40: computes `fig` using `fig=figure('Visible','off')`.
- Lines 41: computes `options.dots` using `options.dots=false`.
- Lines 42: computes `colours` using `colours=(1:4).'`.
- Lines 44: computes `patch_obj` using `patch_obj=findobj(fig,'Type','patch')`.
- Lines 45: computes `line_obj` using `line_obj=findobj(fig,'Type','line')`.
- Lines 46: computes `colour_data` using `colour_data=local_patch_colours(patch_obj)`.
- Lines 52: computes `'setting options.dots` using `'setting options.dots=false must suppress centre marker plotting')`.
- Lines 60: computes `axes_obj` using `axes_obj=findobj(fig,'Type','axes')`.

### Local helper functions

- Line 72: `local_patch_colours()` — `function colour_data=local_patch_colours(patch_obj)`. Preallocate colour data
  - Representative operation: `colour_data=zeros(numel(patch_obj),1)`.
  - Representative operation: `for n=1:numel(patch_obj)`.

## Outputs

- result -regression test result with explanatory messages
- The test draws a tetrahedral spherical Voronoi tessellation, checks the
- patch count and numeric colour mapping, and verifies optional centre dots.

## Implementation structure

- Tests grid_plot() under offscreen graphics. Syntax:
- result=test_dynamic_grid_plot()
- result -regression test result with explanatory messages
- The test draws a tetrahedral spherical Voronoi tessellation, checks the
- patch count and numeric colour mapping, and verifies optional centre dots.
- Announce the test target
- State the grid_plot target of the test
- Force invisible figures during the test
- Build a regular tetrahedral grid on the unit sphere
- Compute the spherical Voronoi tessellation once
- Draw supplied tessera with numeric colours and no centre dots
- Draw with internally generated tessera and default centre dots

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `grid_plot()`, `get()`, `set()`, `onCleanup()`, `xyz()`, `voronoisphere()`, `figure()`, `findobj()`, `local_patch_colours()`, `test_close()`, `test_true()`, `close()`, `isscalar()`, `strcmp()`, `axes_obj()`.
