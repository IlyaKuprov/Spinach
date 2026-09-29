# tests/kernel/test_plotting_helpers_offscreen.m

Source: [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_plotting_helpers_offscreen.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_plotting_helpers_offscreen.m)

## Purpose

Regression test for offscreen execution of Spinach plotting helpers. It exercises the plotting helpers under invisible figures and checks graphics object creation, axis sizes, returned data arrays, and figure helper side effects without relying on image comparison.

## What is tested

- **Offscreen graphics:** figure handles must remain valid while their figures stay invisible. Figure scaling changes the default width and height by factors `1.25` and `0.75` (`1e-12` tolerance); subplots, legends and colour bars must create graphics objects.
- **NMR spectra:** a real `plot_1d` trace retains its input ordinates (`1e-12`) and normal x-axis; a complex trace creates two lines and a legend. `plot_2d` returns the transpose of its spectrum (`1e-12`), frequency axes of the matching lengths, reversed NMR axes, and contour and colour-bar objects. `stack_2d` produces one patch per input column; `plot_3d` produces at least two patch objects and four axes for the volume and projections.
- **MRI and volume output:** image, phantom and k-space modes collectively create three image objects; k-space image data match the real part of the supplied complex data (`1e-12`). Volume rendering exposes surfaces for a small signed 3D test volume.
- **Spatial geometry:** precomputed Voronoi segments match the expected x/y coordinates and NaN separators, whereas an empty tessellation yields `0×0` arrays. The concentration extrusion produces one side patch with 14 vertices and two cap patches with seven vertices each; their connectivity, neutral cap colours and planar polygon areas are checked. The wall area equals the sum over active cells of absolute concentration times cell perimeter (`1e-12` tolerance).

## Inputs and outputs

- **Inputs**: none. The function takes no arguments.
- **Outputs**: `result` — regression test result structure with explanatory messages, accumulated by the `test_true` and `test_close` assertions across all subtests.

## References

- Source file: [tests/kernel/test_plotting_helpers_offscreen.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_plotting_helpers_offscreen.m) in the Spinach repository.
- Tested helpers: `kfigure`, `scale_figure`, `kgrid`, `ktitle`, `kxlabel`, `kylabel`, `kzlabel`, `klegend`, `kcolourbar`, `ksgtitle`, `plot_1d`, `plot_2d`, `stack_2d`, `plot_3d`, `mri_2d_plot`, `volplot`, `mesh_preplot`, `conc_plot`.
