# examples/fitting/fluoroalkanes/anti_difluoroheptane.m

- Signature: `anti_difluoroheptane()`

## Purpose

Fit the 1H NMR spectrum of anti-3,5-difluoroheptane with respect to J-couplings. The fit also includes a 19F spectrum. Methyl groups are ghosted out because they do not influence the signals in question. See our paper for further details: https://doi.org/doi/10.1021/acs.joc.4c00670. Calculation time: hours.

## Model and fitting workflow

- Loads experimental 19F data and two 1H spectral regions, then normalises each spectrum by its maximum.
- Builds a 23-site spin system at a magnet induction of 11.7464, with 12C, 1H, 19F and ghost (`G`) sites. Chemical shifts are specified for the 1H and 19F sites; selected J-couplings are fitted, while couplings to the ghosted methyl sites are fixed at 7.45.
- Uses a 15-parameter initial guess and `fminsearch` to minimise a least-squares error, with `MaxIter` set to 5000 and `MaxFunEvals` to `Inf`. The error is the squared 19F spectral residual plus twice the squared residual for each 1H region.
- Simulates separate 19F and 1H acquisitions without decoupling. The 19F signal receives exponential apodisation; both 1H signals receive Gaussian apodisation. The signals are Fourier transformed with zero filling, reversed, converted to ppm axes and interpolated onto the experimental axes.
- Plots experimental points against simulated spectra for all three regions during optimisation and displays the parameter values and final answer.