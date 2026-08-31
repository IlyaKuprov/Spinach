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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 38-39: Check consistency; implemented by `grumble(f0,df0,f1,df1,x)`.
- Lines 41-43: Adapt to input; implemented by `if isscalar(f0)&&isscalar(df0)&& isscalar(f1)&&isscalar(df1)&&(~isscalar(x))`.
- Lines 48-49: Preallocate the output; implemented by `y=zeros(size(x))`.
- Lines 51-52: Loop over the entries; implemented by `for n=1:numel(x)`.
- Lines 54-55: Get spline coefficients (x^3 -> x^0); implemented by `c=[ 1 2 1 -2`.
- Lines 60-61: Evaluate at the query point; implemented by `y(n)=c(1)*x(n)^3+c(2)*x(n)^2+c(3)*x(n)+c(4)`.

### Control flow inferred from the code

- Line 42: conditional branch on `isscalar(f0)&&isscalar(df0)&&`.
- Line 52: `for` loop over `n=1:numel(x)`.

### Key state/data transformations

- Lines 44: computes `f0` using `f0=f0*ones(size(x)); df0=df0*ones(size(x))`.
- Lines 45: computes `f1` using `f1=f1*ones(size(x)); df1=df1*ones(size(x))`.
- Lines 49: computes `y` using `y=zeros(size(x))`.
- Lines 55: computes `c` using `c=[ 1 2 1 -2`.
- Lines 61: computes `y(n)` using `y(n)=c(1)*x(n)^3+c(2)*x(n)^2+c(3)*x(n)+c(4)`.

### Local helper functions

- Line 68: `grumble()` — `function grumble(f0,df0,f1,df1,x)`.
  - Representative operation: `if (~isnumeric(f0))||(~isreal(f0))||any(~isfinite(f0),'all')|| (~isnumeric(df0))||(~isreal(df0))||any(~isfinite(df0),'all')|| (~isnumeric(f1))||(~isreal(f1))||any(~isfin…`.
  - Representative operation: `(~isnumeric(df0))||(~isreal(df0))||any(~isfinite(df0),'all')|| (~isnumeric(f1))||(~isreal(f1))||any(~isfinite(f1),'all')|| (~isnumeric(df1))||(~isreal(df1))||any(~isfini…`.

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
