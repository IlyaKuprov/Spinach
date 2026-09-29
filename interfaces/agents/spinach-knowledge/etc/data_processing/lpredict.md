# etc/data_processing/lpredict.m

- MATLAB implementation: [etc/data_processing/lpredict.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/data_processing/lpredict.m)

- Signature: `y=lpredict(x,npcoeffs,npredps)`
- [Spinach Wiki page](https://spindynamics.org/wiki/index.php?title=lpredict.m)
- Source authors: malcolm.lidierth@kcl.ac.uk; ilya.kuprov@weizmann.ac.il

## Purpose

Forward linear prediction: extrapolate a real-valued sampled series by fitting a linear predictor and recursively generating future points.

## Call and inputs

`y=lpredict(x,npcoeffs,npredps)`

- `x` — numeric, real column vector with non-zero standard deviation. Its sample spacing and physical units are not changed or inferred; predictions are in the same units as the input values.
- `npcoeffs` — finite real integer predictor order from 2 through `numel(x)` (the source describes this as greater than 1).
- `npredps` — finite positive real integer number of future data points to predict.

## Numerical mechanism

The function stores the mean and standard deviation of `x`, centres and standardises it, and obtains coefficients using MATLAB's `lpc(x,npcoeffs)`. The first predicted point uses the last `npcoeffs` observed samples; following points use a mixture of predictions and remaining observed samples until the recursion is fully prediction-driven. It then restores the original scale and mean.

## Output and limitations

`y` is a column vector containing exactly `npredps` predicted points. The routine returns a single extrapolated sequence, with no uncertainty estimate or validity horizon; predictions therefore depend on the fitted linear predictor remaining useful beyond the observed data. It calls MATLAB's `lpc` function.
