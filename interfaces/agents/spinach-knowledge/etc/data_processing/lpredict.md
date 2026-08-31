# etc/data_processing/lpredict.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/data_processing/lpredict.m`
- Signature: `y=lpredict(x,npcoeffs,npredps)`
- Total lines: 79

## Purpose

Forward linear prediction. Syntax: y=lpredict(x,npcoeffs,npredps)

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 25-26: Check consistency; implemented by `grumble(x,npcoeffs,npredps)`.
- Lines 28-29: Store and subtract the mean; implemented by `mean_x=mean(x); x=x-mean_x`.
- Lines 31-32: Store and scale by stdev; implemented by `std_x=std(x); x=x/std_x`.
- Lines 34-35: Get linear predictor coefficients; implemented by `a=lpc(x,npcoeffs); lpcoeffs=-a(2:end)`.
- Lines 37-38: Pre-allocate output; implemented by `y=zeros(npredps,1)`.
- Lines 40-41: First value only uses x; implemented by `y(1)=lpcoeffs*x(end:-1:(end-npcoeffs+1))`.
- Lines 43-44: Then some x and some y; implemented by `for k=2:min(npcoeffs,npredps)`.
- Lines 49-50: Then keep going using y; implemented by `for k=(npcoeffs+1):npredps`.
- Lines 54-55: Scale and shift back; implemented by `y=std_x*y; y=y+mean_x`.

### Control flow inferred from the code

- Line 44: `for` loop over `k=2:min(npcoeffs,npredps)`.
- Line 50: `for` loop over `k=(npcoeffs+1):npredps`.

### Key state/data transformations

- Lines 29: computes `mean_x` using `mean_x=mean(x); x=x-mean_x`.
- Lines 32: computes `std_x` using `std_x=std(x); x=x/std_x`.
- Lines 35: computes `a` using `a=lpc(x,npcoeffs); lpcoeffs=-a(2:end)`.
- Lines 38: computes `y` using `y=zeros(npredps,1)`.
- Lines 41: computes `y(1)` using `y(1)=lpcoeffs*x(end:-1:(end-npcoeffs+1))`.
- Lines 45-46: computes `y(k)` using `y(k)=lpcoeffs*[y((k-1):-1:1); x(end:-1:(end-npcoeffs+k))]`.

### Local helper functions

- Line 60: `grumble()` — `function grumble(x,npcoeffs,npredps)`.
  - Representative operation: `if (~isnumeric(x))||(~isreal(x))||(~iscolumn(x))||(std(x)==0)`.
  - Representative operation: `error('x must be a real column vector with a non-zero standard deviation.')`.

## Parameters / inputs

- x -input data, a column vector
- npcoeffs -number of predictor coefficients,
- must be greater than 1
- npredps -number of data points to predict

## Outputs

- y -predicted points

## Implementation structure

- Forward linear prediction. Syntax:
- y=lpredict(x,npcoeffs,npredps)
- x -input data, a column vector
- npcoeffs -number of predictor coefficients,
- must be greater than 1
- npredps -number of data points to predict
- y -predicted points
- Check consistency
- Store and subtract the mean
- Store and scale by stdev
- Get linear predictor coefficients
- Pre-allocate output

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `std()`, `lpc()`, `iscolumn()`, `isscalar()`.
