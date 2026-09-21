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
