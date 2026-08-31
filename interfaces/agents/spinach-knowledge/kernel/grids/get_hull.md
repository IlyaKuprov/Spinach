# kernel/grids/get_hull.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/grids/get_hull.m`
- Signature: `[hull,edges]=get_hull(theta_angles,phi_angles)`
- Total lines: 72

## Purpose

Generates a convex hull of a two-angle grid for 2D surface plotting. Syntax: [hull,edges]=get_hull(theta_angles,phi_angles)

## Physical / mathematical content

- Quadrature and geometry utilities. These files generate spherical/SO(3) grids, Voronoi weights, and adaptive integration tools for orientation averaging.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 32-33: Check consistency; implemented by `grumble(theta_angles,phi_angles)`.
- Lines 35-36: Get Cartesian coordinates; implemented by `x=sin(theta_angles).*cos(phi_angles)`.
- Lines 40-41: Get the convex hull; implemented by `hull=convhull(x,y,z)`.
- Lines 43-44: Get the edges; implemented by `edges=unique([hull(:,1) hull(:,2)`.

### Key state/data transformations

- Lines 36: computes `x` using `x=sin(theta_angles).*cos(phi_angles)`.
- Lines 37: computes `y` using `y=sin(theta_angles).*sin(phi_angles)`.
- Lines 38: computes `z` using `z=cos(theta_angles)`.
- Lines 41: computes `hull` using `hull=convhull(x,y,z)`.
- Lines 44: computes `edges` using `edges=unique([hull(:,1) hull(:,2)`.

### Local helper functions

- Line 53: `grumble()` — `function grumble(theta_angles,phi_angles)`.
  - Representative operation: `if (~isnumeric(theta_angles))||(~isreal(theta_angles))|| any(~isfinite(theta_angles))||(size(theta_angles,2)~=1)`.
  - Representative operation: `any(~isfinite(theta_angles))||(size(theta_angles,2)~=1)`.

## Parameters / inputs

- theta_angles -a column vector of theta angles,
- polar coordinates, ISO convention,
- radians
- phi_angles -a column vector of phi angles,
- polar coordinates, ISO convention,
- radians

## Outputs

- hull -a matrix of point indices of
- dimension Nx3, where N is the
- number of triangular facets
- edges -a matrix of point indices of
- dimension Nx2, where N is the
- number of grid edges

## Implementation structure

- Generates a convex hull of a two-angle grid for 2D
- surface plotting. Syntax:
- [hull,edges]=get_hull(theta_angles,phi_angles)
- theta_angles -a column vector of theta angles,
- polar coordinates, ISO convention,
- radians
- phi_angles -a column vector of phi angles,
- hull -a matrix of point indices of
- dimension Nx3, where N is the
- number of triangular facets
- edges -a matrix of point indices of
- dimension Nx2, where N is the

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `convhull()`, `hull()`, `edges()`, `any()`.
