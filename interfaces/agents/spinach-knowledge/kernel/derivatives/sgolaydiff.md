# kernel/derivatives/sgolaydiff.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/derivatives/sgolaydiff.m`
- Signature: `dy=sgolaydiff(y,der_order,npoints,poly_order)`
- Total lines: 111

## Purpose

Savitzky-Golay differentiation of noisy sampled signals by local least-squares polynomial fitting. Syntax: dy=sgolaydiff(y,der_order,npoints,poly_order)

## Physical / mathematical content

- Derivative utilities. These routines compute finite-difference, analytical, or optimisation-oriented derivatives needed for sensitivity analysis, fitting, and optimal control.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 34-35: Check consistency; implemented by `grumble(y,der_order,npoints,poly_order)`.
- Lines 37-38: Count the samples; implemented by `nsamps=size(y,1); dy=0*y`.
- Lines 40-41: Compute sided and centred local fits; implemented by `win_half=(npoints-1)/2; fac=factorial(der_order)`.
- Lines 44-45: Select a local window around the current sample; implemented by `win_left=max(1,min(n-win_half,nsamps-npoints+1))`.
- Lines 48-49: Centre and scale the local sample offsets; implemented by `loc_grid=win_idx(:)-n`.
- Lines 53-54: Build the local Vandermonde matrix; implemented by `V=ones(npoints,poly_order+1)`.
- Lines 59-60: Fit the local polynomial by QR-backed least squares; implemented by `coeff=V\y(win_idx,:)`.
- Lines 62-63: Extract the requested derivative at the expansion point; implemented by `dy(n,:)=fac*coeff(der_order+1,:)/(grid_scale^der_order)`.

### Control flow inferred from the code

- Line 42: `for` loop over `n=1:nsamps`.
- Line 55: `for` loop over `k=1:poly_order`.

### Key state/data transformations

- Lines 38: computes `nsamps` using `nsamps=size(y,1); dy=0*y`.
- Lines 41: computes `win_half` using `win_half=(npoints-1)/2; fac=factorial(der_order)`.
- Lines 45: computes `win_left` using `win_left=max(1,min(n-win_half,nsamps-npoints+1))`.
- Lines 46: computes `win_idx` using `win_idx=win_left:(win_left+npoints-1)`.
- Lines 49: computes `loc_grid` using `loc_grid=win_idx(:)-n`.
- Lines 50: computes `grid_scale` using `grid_scale=max(abs(loc_grid))`.
- Lines 54: computes `V` using `V=ones(npoints,poly_order+1)`.
- Lines 56: computes `V(:,k+1)` using `V(:,k+1)=loc_grid.^k`.
- Lines 60: computes `coeff` using `coeff=V\y(win_idx,:)`.
- Lines 63: computes `dy(n,:)` using `dy(n,:)=fac*coeff(der_order+1,:)/(grid_scale^der_order)`.

### Local helper functions

- Line 70: `grumble()` — `function grumble(y,der_order,npoints,poly_order)`.
  - Representative operation: `if (~isfloat(y))||isempty(y)||issparse(y)||(~ismatrix(y))`.
  - Representative operation: `error('y must be a non-empty dense floating-point matrix.')`.

## Parameters / inputs

- y -N-by-M signal matrix; rows are samples and
- columns are independent signals
- der_order -derivative order; order 0 returns the
- smoothed signal
- npoints -odd number of points in the local least-
- squares window
- poly_order -order of the local polynomial

## Outputs

- dy -N-by-M derivative matrix on a unit-step
- uniform grid
- Note: sgolaydiff(s,1,7,3) is recommended for differentiating
- EPR spectra; use a tight integration tolerance and in-
- crease the number of field/frequency axis points.

## Implementation structure

- Savitzky-Golay differentiation of noisy sampled signals by local
- least-squares polynomial fitting. Syntax:
- dy=sgolaydiff(y,der_order,npoints,poly_order)
- y -N-by-M signal matrix; rows are samples and
- columns are independent signals
- der_order -derivative order; order 0 returns the
- smoothed signal
- npoints -odd number of points in the local least-
- squares window
- poly_order -order of the local polynomial
- dy -N-by-M derivative matrix on a unit-step
- uniform grid

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `factorial()`, `win_idx()`, `coeff()`, `isfloat()`, `issparse()`, `ismatrix()`, `any()`.
