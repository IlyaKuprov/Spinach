# interfaces/comsol/mesh_plot.m

[Canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/comsol/mesh_plot.m) · [Spinach Wiki: mesh_plot.m](https://spindynamics.org/wiki/index.php?title=mesh_plot.m)

## Purpose and return value

`mesh_plot(spin_system,qscale,nodelabels)` draws a prepared two-dimensional mesh in the current plotting axes and returns no output argument. Run `mesh_preplot` first to populate the plotting arrays. The body does not explicitly open a new figure; it calls `hold on`, `box on`, `grid on`, and `axis equal` on the current axes.

## Accepted data and dependencies

`spin_system.mesh` must contain the mesh coordinates `x` and `y`, prepared `plot` data, and (when `qscale` is nonzero) velocity components `u` and `v`. The plotting data used here are `tri_a`/`tri_b`, `rec_a`/`rec_b`, and `vor_a`/`vor_b`. The routine also uses MATLAB graphics (`patch`, `quiver`, `text`) and Spinach's `kxlabel` and `kylabel` axis-label helpers. Its checks confirm that `mesh` and `mesh.plot` exist, but do not validate each required field.

## Drawing, scaling, and guardrails

It draws triangle, rectangle, and Voronoi outlines. The edge-array `edg_a`/`edg_b` patch call is commented out, so those edges are not drawn by this routine. Axis labels read “X position, mm” and “Y position, mm”; coordinates are passed through without unit conversion. If `qscale` is positive, `quiver` receives the mesh coordinates, velocities, and that scale multiplier; zero suppresses velocity arrows. Set `nodelabels` to `1` to label vertices with their row indices or `0` to omit labels. The guards require a non-negative real scalar `qscale` and a real scalar `nodelabels` equal to `0` or `1`. The function changes plotting state and draws graphics; there is no separate success/status return.
