# kernel/plotting/plot_orisel.m

Draws spherical Voronoi polygons colored by a selected row of orientation-selection intensity. The grid has column-vector angles and weights; selection_weights has one column per orientation. Directions use beta and gamma. A full sphere of distinct directions is required; repeated directions in three-angle grids are not supported.

The displayed row is normalized to its maximum (all-zero rows remain zero), with true RGB colors from parula and a [0,1] color scale. Quadrature weights are not applied. Spatial axes and ticks remain visible in a three-dimensional view, and the caller's hold state is restored.
