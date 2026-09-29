# interfaces/comsol/conc_plot.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/comsol/conc_plot.m) · [Spinach Wiki: conc_plot.m](https://spindynamics.org/wiki/index.php?title=conc_plot.m)

## Purpose and call

`conc_plot(spin_system,conc,obs)` draws a vertical bar over each active Voronoi cell, using the supplied concentration value as the bar top and zero as its base. The source describes it as a plotting step to call after `mesh_plot` has drawn the mesh. `obs` is optional.

## Accepted data and geometry

`spin_system.mesh` must contain finite real two-element `zext` data and Voronoi data with `ncells`, `cells`, `vertices`, and `max_cell_size`. `conc` must be a finite real column vector with one value per Voronoi cell. A cell is drawn only when `abs(conc(n)) > 1e-3 * diff(spin_system.mesh.zext)`. Its vertical coordinate is the supplied `conc(n)`; the routine applies no concentration-unit conversion.

When supplied, `obs` must be a finite real matrix with one row per cell and one, two, or three columns. The source treats its columns as phase, then optional amplitude, then optional longitudinal observable. Phase is wrapped with `wrapTo2Pi` and mapped to HSV hue. With one phase column, saturation and value are fixed at 0.75 and 0.50. With phase and amplitude, saturation is amplitude divided by the column maximum (or zero if that maximum is zero), and value is 0.50. With three columns, the same hue and saturation rules apply, while value is the longitudinal column scaled from its minimum-to-maximum range; a constant column maps to value 1. With no `obs`, cells use neutral mid-grey.

## Output and guardrails

There is no returned data structure or status value: the function builds top, bottom, and side faces and draws a flat-coloured patch in the current plotting axes. It checks for the required mesh and Voronoi fields, the shape and finiteness of `zext`, the finite real concentration column and its length, and the observable matrix's finiteness, row count, and column count.
