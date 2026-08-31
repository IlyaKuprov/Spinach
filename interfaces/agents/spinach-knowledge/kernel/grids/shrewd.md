# kernel/grids/shrewd.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/grids/shrewd.m`
- Signature: `weights=shrewd(alphas,betas,gammas,max_rank,max_error)`
- Total lines: 128

## Purpose

Computes SHREWD weights for a given two-or three-angle spherical grid. See the paper by Eden and Levitt for details on now the al- gorithm works: http://dx.doi.org/10.1006/jmre.1998.1427 Syntax: weights=shrewd(alphas,betas,gammas,max_rank,max_error)

## Physical / mathematical content

- Quadrature and geometry utilities. These files generate spherical/SO(3) grids, Voronoi weights, and adaptive integration tools for orientation averaging.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 36-37: Check consistency; implemented by `grumble(alphas,betas,gammas,max_rank,max_error)`.
- Lines 39-40: Decide the grid type; implemented by `if all(alphas==0)`.
- Lines 42-43: Preallocate spherical harmonic matrix; implemented by `H=zeros(lm2lin(max_rank,-max_rank)+1,numel(alphas),'like',1i)`.
- Lines 45-46: Fill spherical harmonic matrix; implemented by `for k=1:numel(alphas)`.
- Lines 55-56: Get the right hand side vector; implemented by `v=max_error*ones(lm2lin(max_rank,-max_rank)+1,1)`.
- Lines 61-62: Preallocate Wigner function matrix; implemented by `H=zeros(lmn2lin(max_rank,-max_rank,-max_rank),numel(alphas))`.
- Lines 64-65: Fill Wigner function matrix; implemented by `for k=1:numel(alphas)`.
- Lines 76-77: Get the right hand side vector; implemented by `v=max_error*ones(lmn2lin(max_rank,-max_rank,-max_rank),1)`.
- Lines 82-83: Compute the weights; implemented by `weights=real(H\v); weights=weights/sum(weights)`.
- Lines 85-86: Run some diagnostics; implemented by `if any(weights==0)`.

### Control flow inferred from the code

- Line 40: conditional branch on `all(alphas==0)`.
- Line 46: `for` loop over `k=1:numel(alphas)`.
- Line 47: `for` loop over `l=0:max_rank`.
- Line 49: `for` loop over `m=l:-1:-l`.
- Line 65: `for` loop over `k=1:numel(alphas)`.
- Line 66: `for` loop over `l=0:max_rank`.
- Line 68: `for` loop over `m=l:-1:-l`.
- Line 69: `for` loop over `n=l:-1:-l`.
- Line 86: conditional branch on `any(weights==0)`.
- Line 89: conditional branch on `any(weights<0)`.

### Key state/data transformations

- Lines 43: computes `H` using `H=zeros(lm2lin(max_rank,-max_rank)+1,numel(alphas),'like',1i)`.
- Lines 48: computes `D` using `D=wigner(l,alphas(k),betas(k),gammas(k))`.
- Lines 50: computes `H(lm2lin(l,m)+1,k)` using `H(lm2lin(l,m)+1,k)=D(l+1,l+m+1)`.
- Lines 56: computes `v` using `v=max_error*ones(lm2lin(max_rank,-max_rank)+1,1)`.
- Lines 57: computes `v(1)` using `v(1)=1-max_error`.
- Lines 70: computes `H(lmn2lin(l,m,n),k)` using `H(lmn2lin(l,m,n),k)=D(l+m+1,l+n+1)`.
- Lines 83: computes `weights` using `weights=real(H\v); weights=weights/sum(weights)`.

### Local helper functions

- Line 96: `grumble()` — `function grumble(alphas,betas,gammas,max_rank,max_error)`.
  - Representative operation: `if (~isnumeric(alphas))||(~isreal(alphas))|| any(~isfinite(alphas))||(size(alphas,2)~=1)`.
  - Representative operation: `any(~isfinite(alphas))||(size(alphas,2)~=1)`.

## Parameters / inputs

- alphas -alpha Euler angles (ZYZ active) of the
- grid, in radians, set to all-zeros for
- two-angle grids
- betas -beta Euler angles (ZYZ active) of the
- grid, in radians
- gammas -gamma Euler angles (ZYZ active) of the
- grid, in radians
- max_rank -maximum spherical rank to take into consi-
- deration when minimizing residuals
- max_error -maximum residual absolute error per spheri-
- cal function

## Outputs

- weights -a vector of grid weights for each
- [alpha beta gamma] point supplied.

## Implementation structure

- Computes SHREWD weights for a given two-or three-angle spherical
- grid. See the paper by Eden and Levitt for details on now the al-
- gorithm works: http://dx.doi.org/10.1006/jmre.1998.1427 Syntax:
- weights=shrewd(alphas,betas,gammas,max_rank,max_error)
- alphas -alpha Euler angles (ZYZ active) of the
- grid, in radians, set to all-zeros for
- two-angle grids
- betas -beta Euler angles (ZYZ active) of the
- grid, in radians
- gammas -gamma Euler angles (ZYZ active) of the
- max_rank -maximum spherical rank to take into consi-
- deration when minimizing residuals

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `all()`, `lm2lin()`, `wigner()`, `alphas()`, `betas()`, `gammas()`, `lmn2lin()`, `any()`.
