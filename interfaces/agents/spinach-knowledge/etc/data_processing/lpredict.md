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
