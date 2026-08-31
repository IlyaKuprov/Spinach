# kernel/plotting/cylgrid.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/cylgrid.m`
- Signature: `cylgrid(zmin,zmax,rmax)`
- Total lines: 99

## Purpose

Draws a cylindrical grid with 10% spacing added around the indicated data extent values. Syntax: cylgrid(zmin,zmax,rmax)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 25-26: Check consistency; implemented by `grumble(zmin,zmax,rmax)`.
- Lines 28-29: Get the extent gaps; implemented by `rgap=0.1*rmax; zgap=0.1*(zmax-zmin)`.
- Lines 31-32: Draw the spokes; implemented by `for n=0:30:330`.
- Lines 48-49: Draw the circles; implemented by `phi_grid=linspace(0,2*pi,360)`.
- Lines 55-58: Set axis extents; implemented by `axis([(-rmax-2*rgap) (rmax+2*rgap) (-rmax-2*rgap) (rmax+2*rgap) (+zmin-2*zgap) (zmax+2*zgap)])`.
- Lines 60-61: Draw vertical tick labels; implemented by `tick_vals=linspace(zmin-zgap,zmax+zgap,7)`.
- Lines 69-71: Switch on the perspective, kill default axes; implemented by `set(gca,'Projection','perspective','Box','off', 'XTick',[],'YTick',[],'ZTick',[],'visible','off')`.

### Control flow inferred from the code

- Line 32: `for` loop over `n=0:30:330`.
- Line 50: `for` loop over `r=linspace(0,rmax+rgap,7)`.
- Line 62: `for` loop over `n=1:numel(tick_vals)`.

### Key state/data transformations

- Lines 29: computes `rgap` using `rgap=0.1*rmax; zgap=0.1*(zmax-zmin)`.
- Lines 49: computes `phi_grid` using `phi_grid=linspace(0,2*pi,360)`.
- Lines 61: computes `tick_vals` using `tick_vals=linspace(zmin-zgap,zmax+zgap,7)`.

### Local helper functions

- Line 77: `grumble()` — `function grumble(zmin,zmax,rmax)`.
  - Representative operation: `if (~isnumeric(zmin))||(~isreal(zmin))|| (numel(zmin)~=1)||(~isfinite(zmin))`.
  - Representative operation: `(numel(zmin)~=1)||(~isfinite(zmin))`.

## Parameters / inputs

- zmin -lower bound on the Z axis
- zmax -upper bound on the Z axis
- rmax -upper bound on the radius

## Outputs

- this function updates the current figure

## Implementation structure

- Draws a cylindrical grid with 10% spacing added around
- the indicated data extent values. Syntax:
- cylgrid(zmin,zmax,rmax)
- zmin -lower bound on the Z axis
- zmax -upper bound on the Z axis
- rmax -upper bound on the radius
- this function updates the current figure
- Check consistency
- Get the extent gaps
- Draw the spokes
- Draw the circles
- Set axis extents

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `line()`, `cosd()`, `sind()`, `text()`, `num2str()`, `tick_vals()`, `set()`.
