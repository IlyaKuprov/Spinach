# kernel/utilities/jacobianest.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/jacobianest.m`
- Signature: `[jac,err] = jacobianest(fun,x0)`
- Total lines: 177

## Purpose

Estimate of the Jacobian matrix of a vector valued function of n variables. Syntax: [jac,err] = jacobianest(fun,x0)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `swapelement()`, `rombextrap()`, `vec2mat()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 31-32: Check consistency; implemented by `grumble(fun,x0)`.
- Lines 34-35: get the length of x0 for the size of jac; implemented by `nx = numel(x0)`.
- Lines 37-38: Finite difference gridding; implemented by `MaxStep = 100`.
- Lines 43-44: was a string supplied?; implemented by `if ischar(fun)`.
- Lines 48-49: get fun at the center point; implemented by `f0 = fun(x0)`.
- Lines 53-54: empty begets empty; implemented by `jac = zeros(0,nx)`.
- Lines 59-60: total number of derivatives we will need to take; implemented by `jac = zeros(n,nx)`.
- Lines 70-72: evaluate at each step, centered around x0_i difference to give a second order estimate; implemented by `fdel = zeros(n,nsteps)`.
- Lines 79-81: these are pure second order estimates of the first derivative, for each trial delta.; implemented by `derest = fdel.*repmat(0.5 ./ delta,n,1)`.
- Lines 83-90: The error term on these estimates has a second order component, but also some 4th and 6th order terms in it. Use Romberg exrapolation to improve the estimates to 6th order, as well as to provide the error estimate.; implemented by `for j = 1:n`.
- Lines 88-90: loop here, as rombextrap coupled with the trimming will get complicated otherwise.; implemented by `for j = 1:n`.
- Lines 93-94: trim off 3 estimates at each end of the scale; implemented by `nest = length(der_romb)`.
- Lines 101-102: now pick the estimate with the lowest predicted error; implemented by `[err(j,i),ind] = min(errest)`.

### Control flow inferred from the code

- Line 44: conditional branch on `ischar(fun)`.
- Line 52: conditional branch on `n==0`.
- Line 62: `for` loop over `i = 1:nx`.
- Line 64: conditional branch on `x0_i ~= 0`.
- Line 73: `for` loop over `j = 1:nsteps`.
- Line 90: `for` loop over `j = 1:n`.

### Key state/data transformations

- Lines 35: computes `nx` using `nx = numel(x0)`.
- Lines 38: computes `MaxStep` using `MaxStep = 100`.
- Lines 39: computes `StepRatio` using `StepRatio = 2.0000001`.
- Lines 40: computes `relativedelta` using `relativedelta = MaxStep*StepRatio .^(0:-1:-25)`.
- Lines 41: computes `nsteps` using `nsteps = length(relativedelta)`.
- Lines 45: computes `fun` using `fun = str2func(fun)`.
- Lines 49: computes `f0` using `f0 = fun(x0)`.
- Lines 51: computes `n` using `n = length(f0)`.
- Lines 54: computes `jac` using `jac = zeros(0,nx)`.
- Lines 55: computes `err` using `err = jac`.
- Lines 63: computes `x0_i` using `x0_i = x0(i)`.
- Lines 65: computes `delta` using `delta = x0_i*relativedelta`.
- Lines 72: computes `fdel` using `fdel = zeros(n,nsteps)`.
- Lines 74: computes `fdif` using `fdif = fun(swapelement(x0,i,x0_i + delta(j))) - fun(swapelement(x0,i,x0_i - delta(j)))`.
- Lines 76: computes `fdel(:,j)` using `fdel(:,j) = fdif(:)`.
- Lines 81: computes `derest` using `derest = fdel.*repmat(0.5 ./ delta,n,1)`.
- Lines 91: computes `[der_romb,errest]` using `[der_romb,errest] = rombextrap(StepRatio,derest(j,:),[2 4])`.
- Lines 94: computes `nest` using `nest = length(der_romb)`.

### Local helper functions

- Line 111: `swapelement()` — `function vec = swapelement(vec,ind,val)`. Romberg extrapolation StepRatio -Ratio decrease in step
  - Representative operation: `vec(ind) = val`.
- Line 125: `rombextrap()` — `function [der_romb,errest] = rombextrap(StepRatio,der_init,rombexpon)`. Inverse step ratio
  - Representative operation: `srinv = 1/StepRatio`.
  - Representative operation: `nexpon = length(rombexpon)`.
- Line 158: `vec2mat()` — `function mat = vec2mat(vec,n,m)`. Consistency enforcement
  - Representative operation: `[i,j] = ndgrid(1:n,0:m-1)`.
  - Representative operation: `ind = i+j; mat = vec(ind)`.
- Line 165: `grumble()` — `function grumble(fun,x0)`. To prove that you are not a robot, injure a human being, or, through inaction, allow a human being
  - Representative operation: `if ~isa(fun,'function_handle')`.
  - Representative operation: `error('fun must be a function handle.')`.

## Parameters / inputs

- fun -(vector valued) analytical function to differentiate.
- fun must be a function of the vector or array x0.
- x0 -vector location at which to differentiate fun
- If x0 is an nxm array, then fun is assumed to be
- a function of n*m variables.

## Outputs

- jac -array of first partial derivatives of fun.
- Assuming that x0 is a vector of length p
- and fun returns a vector of length n, then
- jac will be an array of size (n,p)
- err -vector of error estimates corresponding to
- each partial derivative in jac.
- John D'Errico

## Implementation structure

- Estimate of the Jacobian matrix of a vector valued
- function of n variables. Syntax:
- [jac,err] = jacobianest(fun,x0)
- fun -(vector valued) analytical function to differentiate.
- fun must be a function of the vector or array x0.
- x0 -vector location at which to differentiate fun
- If x0 is an nxm array, then fun is assumed to be
- a function of n*m variables.
- jac -array of first partial derivatives of fun.
- Assuming that x0 is a vector of length p
- and fun returns a vector of length n, then
- jac will be an array of size (n,p)

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ischar()`, `str2func()`, `fun()`, `swapelement()`, `delta()`, `fdel()`, `fdif()`, `rombextrap()`, `derest()`, `der_romb()`, `tags()`, `errest()`, `err()`, `jac()`, `vec()`.
