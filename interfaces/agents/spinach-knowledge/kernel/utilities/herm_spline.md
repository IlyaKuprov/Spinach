# kernel/utilities/herm_spline.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/herm_spline.m`
- Signature: `y=herm_spline(f0,df0,f1,df1,x)`
- Total lines: 115

## Purpose

Cubic Hermite spline on [0,1] interval from values and deriva- tives at the interval edges. Syntax: y=herm_spline(f0,df0,f1,df1,x)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- f0 -function value(s) at the left edge, a real
- scalar or array
- df0 -function derivative(s) at the left edge,
- a real scalar or array
- f1 -function value(s) at the right edge, a real
- scalar or array
- df1 -function derivative(s) at the right edge,
- a real scalar or array
- x -query point(s) inside [0,1] interval,
- a real scalar or array
- Function values and derivatives can be scalars (in which case
- the same spline is evaluated at all query points) or arrays of
- the same size as x, in which case multiple splines are evalua-
- ted at their corresponding query points.

## Outputs

- y -the value of the spline(s) at the query point(s)

## Implementation structure

- Cubic Hermite spline on [0,1] interval from values and deriva-
- tives at the interval edges. Syntax:
- y=herm_spline(f0,df0,f1,df1,x)
- f0 -function value(s) at the left edge, a real
- scalar or array
- df0 -function derivative(s) at the left edge,
- a real scalar or array
- f1 -function value(s) at the right edge, a real
- df1 -function derivative(s) at the right edge,
- x -query point(s) inside [0,1] interval,
- Function values and derivatives can be scalars (in which case
- the same spline is evaluated at all query points) or arrays of

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isscalar()`, `df0()`, `df1()`, `any()`, `all()`.
