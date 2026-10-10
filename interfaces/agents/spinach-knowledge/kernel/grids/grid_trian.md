# kernel/grids/grid_trian.m

- Signature: `[alps,bets,gams,whts,vorn]=grid_trian(type,n)`

## Behaviour

Generates a spherical triangular quadrature grid. `type` is a character string selecting `'asg'`, `'sophe'`, or `'stoll'`; `n` is a positive real integer subdivision parameter, not a returned point count. All three outputs `alps`, `bets`, and `gams` are matching column vectors with one entry per grid point, in radians; `alps` is zero throughout because these are two-angle grids.

The `'asg'` branch maps an integer triangular lattice onto the sphere and extends it by coordinate reflections. `'sophe'` samples a triangular beta-gamma lattice in an octant, reflects it over the sphere, and explicitly adds the two poles. `'stoll'` combines three SOPHE-derived octant constructions, reflects the result, and adds the poles and four equatorial points. The number and arrangement of points depend on the selected construction and `n`.

When more than three outputs are requested, the function computes spherical Voronoi polygons and solid-angle cell areas; `whts` is the area vector divided by `4*pi`, so its entries are normalised solid-angle weights. `vorn` contains one spherical Voronoi polygon per grid point. This calculation is skipped for calls requesting only the three angle vectors. With no output arguments the function also computes the tessellation and plots the grid.

The source identifies Appendix A.6 of the cited paper as the grid reference: [DOI](http://dx.doi.org/10.1016/j.jmr.2014.05.009).

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/grids/grid_trian.m) · [Wiki](https://spindynamics.org/wiki/index.php?title=grid_trian.m)