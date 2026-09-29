# tests/kernel/test_dynamic_grid_plot.m

Source: [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_grid_plot.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_grid_plot.m)

## Purpose

Regression test for `grid_plot()` under offscreen graphics. It verifies that `grid_plot()` draws one patch per spherical Voronoi cell, passes numeric colour data through to the patch objects, honours the `options.dots` centre-dot setting, and leaves the axes in square plotting mode.

## Behaviour

- Announces the test target with `fprintf('TESTING: Spherical grid plotting\n')` and initialises the result via `new_test_result()` for the `kernel/dynamic_grid_plot` target, described as "Offscreen spherical grid plotting".
- Saves the current root default figure visibility, sets `groot` `defaultFigureVisible` to `'off'`, and registers an `onCleanup` call that restores the original visibility when the test finishes.
- Builds a regular tetrahedral grid on the unit sphere from the vertex matrix `[1 1 -1 -1;1 -1 1 -1;1 -1 -1 1]`, normalising each column by its Euclidean norm, and extracts `x`, `y`, `z` column vectors.
- Computes the spherical Voronoi tessellation once with `voronoisphere(xyz)`.
- First drawing pass: creates an invisible figure, sets `options.dots=false`, defines `colours=(1:4).'`, and calls `grid_plot(x,y,z,vorn,colours,options)`. It then:
  - Checks with `test_close` that the number of patch objects found by `findobj(fig,'Type','patch')` equals `numel(vorn)` with tolerances `0,0`, on the grounds that a tetrahedral Voronoi grid must render one patch per tessellation cell.
  - Extracts scalar `CData(1)` from each patch via the local helper `local_patch_colours` and checks with `test_close` that `sort(colour_data)` matches `colours` with tolerances `1e-15,1e-15`, so numeric colour data is passed through cell by cell.
  - Checks with `test_true` that no line objects exist (`findobj(fig,'Type','line')` is empty), confirming `options.dots=false` suppresses centre marker plotting.
- Second drawing pass: creates another invisible figure and calls `grid_plot(x,y,z)` with tessera, colours, and options omitted. It then:
  - Checks with `test_close` that the generated patch count equals `numel(vorn)` with tolerances `0,0`, i.e. the same cell count is produced when tessera are omitted.
  - Checks with `test_true` that exactly one line object exists, confirming the default options draw one line object containing centre dots.
  - Checks with `test_true` that the first axes object has `PlotBoxAspectRatioMode` equal to `'manual'`, confirming square plotting mode.
- Closes each figure after its checks.

## Inputs and outputs

- Syntax: `result=test_dynamic_grid_plot()`
- Inputs: none.
- Outputs: `result` — regression test result structure with explanatory messages, produced by `new_test_result()` and updated by `test_close` and `test_true`.

## References

- Tested function: `grid_plot()`
- Tessellation helper: `voronoisphere()`
- Test utilities: `new_test_result()`, `test_close()`, `test_true()`
