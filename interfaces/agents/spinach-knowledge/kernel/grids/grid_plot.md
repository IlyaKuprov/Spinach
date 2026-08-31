# kernel/grids/grid_plot.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/grids/grid_plot.m`
- Signature: `grid_plot(x,y,z,vorn,c,options)`
- Total lines: 126

## Purpose

Spherical quadrature grid plotter. Takes a cloud of points on a sphere and plots its Voronoi tessellation. Syntax: grid_plot(x,y,z,vorn,c,options)

## Physical / mathematical content

- Quadrature and geometry utilities. These files generate spherical/SO(3) grids, Voronoi weights, and adaptive integration tools for orientation averaging.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `size()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 32-33: Run Voronoi tessellation on a sphere; implemented by `if (~exist('vorn','var'))||isempty(vorn)`.
- Lines 37-38: Default tessellation face colour is white; implemented by `if (~exist('c','var'))||isempty(c), c='white'; end`.
- Lines 40-42: Default is to draw centre dots; implemented by `if (~exist('options','var'))|| (~isfield(options,'dots'))`.
- Lines 46-47: Check consistency; implemented by `grumble(x,y,z,vorn)`.
- Lines 49-50: Plot the points; implemented by `if options.dots`.
- Lines 56-57: Plot the tessera; implemented by `for k=1:numel(vorn)`.
- Lines 60-62: Colour by the name; implemented by `patch(vorn{k}(1,:),vorn{k}(2,:), vorn{k}(3,:),c,'FaceAlpha',1)`.
- Lines 66-68: Colour by the numbers; implemented by `patch(vorn{k}(1,:),vorn{k}(2,:), vorn{k}(3,:),c(k),'FaceAlpha',1)`.
- Lines 73-74: Residual cosmetics; implemented by `axis square; box on; campos([0 0 10])`.

### Control flow inferred from the code

- Line 33: conditional branch on `(~exist('vorn','var'))||isempty(vorn)`.
- Line 38: conditional branch on `(~exist('c','var'))||isempty(c), c='white'; end`.
- Line 41: conditional branch on `(~exist('options','var'))||`.
- Line 50: conditional branch on `options.dots`.
- Line 57: `for` loop over `k=1:numel(vorn)`.
- Line 58: conditional branch on `ischar(c)`.

### Key state/data transformations

- Lines 34: computes `[~,~,vorn]` using `[~,~,vorn]=voronoisphere([x'; y'; z'])`.
- Lines 43: computes `options.dots` using `options.dots=true`.

### Local helper functions

- Line 80: `grumble()` — `function grumble(x,y,z,vorn)`.
  - Representative operation: `if (~isnumeric(x))||(~isreal(x))||any(~isfinite(x),'all')|| (~isnumeric(y))||(~isreal(y))||any(~isfinite(y),'all')|| (~isnumeric(z))||(~isreal(z))||any(~isfinite(z),'all…`.
  - Representative operation: `(~isnumeric(y))||(~isreal(y))||any(~isfinite(y),'all')|| (~isnumeric(z))||(~isreal(z))||any(~isfinite(z),'all')`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `voronoisphere()`, `isfield()`, `grumble()`, `plot3()`, `xlim()`, `ylim()`, `zlim()`, `ischar()`, `patch()`, `campos()`, `xticks()`, `yticks()`, `zticks()`, `any()`, `iscolumn()`.
