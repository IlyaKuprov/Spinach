# etc/data_processing/lpredict.m

- Signature: `y=lpredict(x,npcoeffs,npredps)`

## Purpose

Forward linear prediction. Syntax: y=lpredict(x,npcoeffs,npredps)

## Physical / mathematical content

## Numerical / algorithmic content

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
