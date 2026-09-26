# kernel/grids/grid_plot.m

- Signature: `grid_plot(x,y,z,vorn,c,options)`

## Purpose

Spherical quadrature grid plotter. Takes a cloud of points on a sphere and plots its Voronoi tessellation. Syntax: grid_plot(x,y,z,vorn,c,options)

## Physical / mathematical content

- Quadrature and geometry utilities. These files generate spherical/SO(3) grids, Voronoi weights, and adaptive integration tools for orientation averaging.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Parameters / inputs

- x,y,z -column vectors containing Cartesian
- coordinates of grid points
- c -values to be mapped into the colour
- of each tessellation face, white if
- this input is left empty
- vorn -Voronoi tessellation; if this is not
- provided, it will be computed
- options.dots -the default (true) puts black
- dots at centres of tessellati-
- on faces

## Outputs

- this function plots a figure

## Implementation structure

- Spherical quadrature grid plotter. Takes a cloud of points
- on a sphere and plots its Voronoi tessellation. Syntax:
- grid_plot(x,y,z,vorn,c,options)
- x,y,z -column vectors containing Cartesian
- coordinates of grid points
- c -values to be mapped into the colour
- of each tessellation face, white if
- this input is left empty
- vorn -Voronoi tessellation; if this is not
- provided, it will be computed
- options.dots -the default (true) puts black
- dots at centres of tessellati-
