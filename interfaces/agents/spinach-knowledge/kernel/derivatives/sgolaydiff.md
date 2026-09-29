# kernel/derivatives/sgolaydiff.m

[Direct MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/derivatives/sgolaydiff.m) · [Spinach Wiki documentation](https://spindynamics.org/wiki/index.php?title=sgolaydiff.m)

## Purpose and signature

`dy=sgolaydiff(y,der_order,npoints,poly_order)` smooths or differentiates sampled signals by local polynomial least-squares fits.

## Inputs and output

- `y`: finite, non-empty, dense floating-point `N`-by-`M` matrix. Rows are samples and columns are independent signals.
- `der_order`: non-negative integer derivative order, no greater than `poly_order`; order zero returns the fitted smoothed values.
- `npoints`: odd integer window width, at least 3 and no greater than `N`.
- `poly_order`: non-negative integer strictly less than `npoints`.
- `dy`: `N`-by-`M` output, with one requested derivative value per input sample and signal.

There is no abscissa or sample-spacing argument: derivatives are with respect to a uniform unit-step sample coordinate. To express derivatives against an axis with spacing `h`, scale the returned order-`der_order` derivative by `h^(-der_order)`.

## Local fits and boundaries

At each row, the routine selects a window of `npoints` samples. Interior windows are centred on the current sample; near either edge the window is shifted to remain inside the data, producing one-sided fits at the boundaries. It centres the integer sample offsets at the current row and scales them by their largest absolute offset, builds a Vandermonde matrix through degree `poly_order`, and solves `V\y(window,:)` by MATLAB least squares. The requested coefficient is multiplied by `factorial(der_order)` and divided by the offset scale to the same derivative order.

The guard also requires finite input data, at least three input rows, a valid odd window, and the stated order/window relationships.

## Existing usage note

`sgolaydiff(s,1,7,3)` is recommended for differentiating EPR spectra; use a tight integration tolerance and increase the number of field/frequency axis points.
