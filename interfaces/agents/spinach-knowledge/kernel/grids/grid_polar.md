# kernel/grids/grid_polar.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/grids/grid_polar.m`
- Signature: `[phi,r,L]=grid_polar(ncircles,rmax)`
- Total lines: 98

## Purpose

Generates a balanced polar grid in which the density of points does not increase towards the centre. Syntax: [phi,r,L]=grid_polar(ncircles,rmax)

## Physical / mathematical content

- Quadrature and geometry utilities. These files generate spherical/SO(3) grids, Voronoi weights, and adaptive integration tools for orientation averaging.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 33-34: Check consistency; implemented by `grumble(ncircles,rmax)`.
- Lines 36-37: Get the radius grid; implemented by `radii=linspace(0,rmax,ncircles)`.
- Lines 39-40: Preallocate output arrays; implemented by `phi=zeros(3*(ncircles-1)^2+2*(ncircles-1)+1,1)`.
- Lines 43-44: Build grid circles; implemented by `for n=1:(numel(radii)-1)`.
- Lines 46-47: Generate the circle; implemented by `phi_current=linspace(0,2*pi,6*n)'; phi_current(end)=[]`.
- Lines 50-51: Update the output array; implemented by `phi((3*(n-1)^2+2*(n-1)+2):(3*n^2+2*n+1))=phi_current`.
- Lines 56-57: Build the Laplacian; implemented by `if nargout>2`.
- Lines 59-60: Make the grid cartesian; implemented by `x=r.*cos(phi); y=r.*sin(phi)`.
- Lines 62-63: Run the triangulation; implemented by `tri=delaunay(x,y)`.
- Lines 65-66: Populate the Laplacian; implemented by `L=zeros(numel(r),numel(r))`.
- Lines 73-74: Do the miscellaneous cosmetics; implemented by `L=L+transpose(L); L=L-diag(sum(L,1)); L=-L/norm(L,2); L=sparse(L)`.

### Control flow inferred from the code

- Line 44: `for` loop over `n=1:(numel(radii)-1)`.
- Line 57: conditional branch on `nargout>2`.
- Line 67: `for` loop over `n=1:size(tri,1)`.

### Key state/data transformations

- Lines 37: computes `radii` using `radii=linspace(0,rmax,ncircles)`.
- Lines 40: computes `phi` using `phi=zeros(3*(ncircles-1)^2+2*(ncircles-1)+1,1)`.
- Lines 41: computes `r` using `r=zeros(3*(ncircles-1)^2+2*(ncircles-1)+1,1)`.
- Lines 47: computes `phi_current` using `phi_current=linspace(0,2*pi,6*n)'; phi_current(end)=[]`.
- Lines 48: computes `r_current` using `r_current=radii(n+1)*ones(6*n,1); r_current(end)=[]`.
- Lines 51: computes `phi((3*(n-1)^2+2*(n-1)+2):(3*n^2+2*n+1))` using `phi((3*(n-1)^2+2*(n-1)+2):(3*n^2+2*n+1))=phi_current`.
- Lines 52: computes `r((3*(n-1)^2+2*(n-1)+2):(3*n^2+2*n+1))` using `r((3*(n-1)^2+2*(n-1)+2):(3*n^2+2*n+1))=r_current`.
- Lines 60: computes `x` using `x=r.*cos(phi); y=r.*sin(phi)`.
- Lines 63: computes `tri` using `tri=delaunay(x,y)`.
- Lines 66: computes `L` using `L=zeros(numel(r),numel(r))`.
- Lines 68: computes `L(tri(n,1),tri(n,2))` using `L(tri(n,1),tri(n,2))=1/norm([x(tri(n,1))-x(tri(n,2)) y(tri(n,1))-y(tri(n,2))],2)^2`.
- Lines 69: computes `L(tri(n,1),tri(n,3))` using `L(tri(n,1),tri(n,3))=1/norm([x(tri(n,1))-x(tri(n,3)) y(tri(n,1))-y(tri(n,3))],2)^2`.
- Lines 70: computes `L(tri(n,2),tri(n,3))` using `L(tri(n,2),tri(n,3))=1/norm([x(tri(n,2))-x(tri(n,3)) y(tri(n,2))-y(tri(n,3))],2)^2`.

### Local helper functions

- Line 81: `grumble()` — `function grumble(ncircles,rmax)`. Буря мглою месит жижу,
  - Representative operation: `if (~isnumeric(ncircles))||(~isreal(ncircles))|| (~isscalar(ncircles))||(ncircles<2)||mod(ncircles,1)`.
  - Representative operation: `(~isscalar(ncircles))||(ncircles<2)||mod(ncircles,1)`.

## Parameters / inputs

- ncircles -number of radial circles
- in the grid, an integer
- rmax -maximum radius that the
- grid must reach

## Outputs

- phi -a column vector of polar
- phi angles, radians
- rmax -a coluim vector of radii
- L -sparse Laplacian operator
- acting on functions defined
- as vector of values in the
- same order as the grid

## Implementation structure

- Generates a balanced polar grid in which the density of
- points does not increase towards the centre. Syntax:
- [phi,r,L]=grid_polar(ncircles,rmax)
- ncircles -number of radial circles
- in the grid, an integer
- rmax -maximum radius that the
- grid must reach
- phi -a column vector of polar
- phi angles, radians
- rmax -a coluim vector of radii
- L -sparse Laplacian operator
- acting on functions defined

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `phi_current()`, `radii()`, `r_current()`, `phi()`, `delaunay()`, `tri()`, `transpose()`, `isscalar()`.
