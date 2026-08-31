# kernel/grids/grid_igloo.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/grids/grid_igloo.m`
- Signature: `[alps,bets,gams,whts,vorn]=grid_igloo(n_long)`
- Total lines: 104

## Purpose

Igloo grid, as per Appendix A2 of

## Physical / mathematical content

- Quadrature and geometry utilities. These files generate spherical/SO(3) grids, Voronoi weights, and adaptive integration tools for orientation averaging.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 36-37: Check consistency; implemented by `grumble(n_long)`.
- Lines 39-40: Betas equispaced by longitude; implemented by `bets=linspace(0,pi,n_long)'`.
- Lines 42-43: Preallocate three-angle blocks; implemented by `alps_bets_gams=cell(n_long,3)`.
- Lines 45-46: Loop over longitudes; implemented by `for k=1:n_long`.
- Lines 48-49: Compute the number of lattitudes; implemented by `n_latt=floor(2*(n_long-1)*sin(bets(k))+0.5)`.
- Lines 52-53: Get the gamma angles of lattitudes; implemented by `gams=linspace(0,2*pi,n_latt+1)'`.
- Lines 56-57: Fill the angle block; implemented by `alps_bets_gams{k,1}=zeros(size(gams))`.
- Lines 63-64: Concatenate angle blocks; implemented by `alps_bets_gams=cell2mat(alps_bets_gams)`.
- Lines 66-67: Extract individual angles; implemented by `alps=alps_bets_gams(:,1)`.
- Lines 71-72: Weights are expensive; implemented by `if (nargout>3)||(nargout==0)`.
- Lines 74-75: Get point coordinates; implemented by `x=sin(bets).*cos(gams)`.
- Lines 79-80: Run Voronoi tessellation; implemented by `[~,~,vorn,whts]=voronoisphere([x'; y'; z'])`.
- Lines 82-83: Get weights; implemented by `whts=whts/(4*pi)`.
- Lines 87-88: If no output requested, plot the grid; implemented by `if nargout==0, grid_plot(x,y,z,vorn); end`.

### Control flow inferred from the code

- Line 46: `for` loop over `k=1:n_long`.
- Line 72: conditional branch on `(nargout>3)||(nargout==0)`.
- Line 88: conditional branch on `nargout==0, grid_plot(x,y,z,vorn); end`.

### Key state/data transformations

- Lines 40: computes `bets` using `bets=linspace(0,pi,n_long)'`.
- Lines 43: computes `alps_bets_gams` using `alps_bets_gams=cell(n_long,3)`.
- Lines 49: computes `n_latt` using `n_latt=floor(2*(n_long-1)*sin(bets(k))+0.5)`.
- Lines 53: computes `gams` using `gams=linspace(0,2*pi,n_latt+1)'`.
- Lines 54: computes `gams(end)` using `gams(end)=[]`.
- Lines 57: computes `alps_bets_gams{k,1}` using `alps_bets_gams{k,1}=zeros(size(gams))`.
- Lines 58: computes `alps_bets_gams{k,2}` using `alps_bets_gams{k,2}=bets(k)*ones(size(gams))`.
- Lines 59: computes `alps_bets_gams{k,3}` using `alps_bets_gams{k,3}=gams`.
- Lines 67: computes `alps` using `alps=alps_bets_gams(:,1)`.
- Lines 75: computes `x` using `x=sin(bets).*cos(gams)`.
- Lines 76: computes `y` using `y=sin(bets).*sin(gams)`.
- Lines 77: computes `z` using `z=cos(bets)`.
- Lines 80: computes `[~,~,vorn,whts]` using `[~,~,vorn,whts]=voronoisphere([x'; y'; z'])`.
- Lines 83: computes `whts` using `whts=whts/(4*pi)`.

### Local helper functions

- Line 93: `grumble()` — `function grumble(n_long)`. Разговор в автобусе: -Он не хочет, говорит что у него иссяк запал.
  - Representative operation: `if (~isnumeric(n_long))||(~isreal(n_long))|| (~isscalar(n_long))||(n_long<1)||mod(n_long,1)`.
  - Representative operation: `(~isscalar(n_long))||(n_long<1)||mod(n_long,1)`.

## Syntax

```matlab
[alps,bets,gams,whts,vorn]=grid_igloo(n_long)
```

## Parameters / inputs

- n_long -number of longitudes in the grid

## Outputs

- alps -alpha Euler angles of the grid (radians),
- zeros because these are two-angle grids
- bets -beta Euler angles of the grid (radians)
- gams -gamma Euler angles of the grid (radians)
- whts -Voronoi tessellation body angle weights
- vorn -a cell array of matrices containing the
- coordinates of the vertices of the Voro-
- noi polyhedra
- If no outputs are requested, a schematic is drawn.

## Implementation structure

- Igloo grid, as per Appendix A2 of
- [alps,bets,gams,whts,vorn]=grid_igloo(n_long)
- n_long -number of longitudes in the grid
- alps -alpha Euler angles of the grid (radians),
- zeros because these are two-angle grids
- bets -beta Euler angles of the grid (radians)
- gams -gamma Euler angles of the grid (radians)
- whts -Voronoi tessellation body angle weights
- vorn -a cell array of matrices containing the
- coordinates of the vertices of the Voro-
- noi polyhedra
- If no outputs are requested, a schematic is drawn.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `bets()`, `gams()`, `cell2mat()`, `alps_bets_gams()`, `voronoisphere()`, `grid_plot()`, `isscalar()`.
