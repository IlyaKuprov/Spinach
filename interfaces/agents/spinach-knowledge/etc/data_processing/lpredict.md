# etc/data_processing/lpredict.m

- Signature: `y=lpredict(x,npcoeffs,npredps)`
- [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=lpredict.m)

## Purpose and method

Predict future samples of a real-valued time series by linear prediction. The routine centres and standardises the input, estimates autoregressive coefficients with MATLAB’s `lpc`, then recursively generates the requested samples. It rescales and recentres the predictions before returning them.

## Inputs

- `x` — real column vector with non-zero standard deviation.
- `npcoeffs` — integer predictor order from 2 through `numel(x)`.
- `npredps` — positive integer number of samples to predict.

## Output

- `y` — column vector of `npredps` predicted samples.
