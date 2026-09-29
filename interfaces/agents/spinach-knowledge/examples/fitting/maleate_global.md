# examples/fitting/maleate_global.m

- MATLAB implementation: [examples/fitting/maleate_global.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fitting/maleate_global.m)

- Signature: `maleate_global()`

## Purpose

Jointly fits the proton and carbon spectra of a slightly asymmetric maleate diester. The source estimates hours for this example; it does not report a fitted result.

## Data and fit vector

The function expects `maleate_1h.mat` and `maleate_13c.mat` in the MATLAB working directory. Each file supplies `axis_hz` and `data_expt`; the script retains the axes and normalises each spectrum by its own maximum. A shared 12-element vector sets `c_shift_1`, `c_shift_2`, `h_shift_1`, `h_shift_2`, then `j_cc`, far and near `j_ch`, and `j_hh` entries, followed by proton/carbon signal scales `a_h`, `a_c` and apodisation settings `lw_h`, `lw_c`. The initial guess is `[-0.0560, 0.0358, 0.0005, -0.0013, 71.6099, -1.1646, 166.5465, 11.9548, 0.0385, 0.0344, 10.1085, 6.2780]`; it is not a fitted result. The source assigns these values to shift/coupling/apodisation fields but does not state their units.

## Spin and acquisition model

The trial vector is used for both nuclei in a four-spin `1H/13C/13C/1H` liquid-state model. The source sets the field to `2*pi*500.101412e6/spin('1H')` and uses a Zeeman-Hilbert basis without approximation. Separate `@acquire` simulations use zero offset, no decoupling, 600 sweep, 1024 acquired points, 4096-point zero fill, and Hz axes for proton and carbon. Their FIDs receive independent exponential apodisation and division by 10; Fourier-transformed real spectra are scaled by `a_h` and `a_c` respectively. The theoretical axes come from `sweep2ticks`, and each simulated spectrum is shape-preserving-cubic interpolated with `interp1(...,'pchip')` onto its experimental axis before plotting and residual evaluation.

## Objective and what the entry point exposes

`maleate_global()` runs `fminsearch` from the initial vector, with iterative display, a 5000-iteration cap, and an unlimited function-evaluation setting. The local objective is the sum of squared residual norms against the separately normalised proton and carbon data. It draws a two-panel experimental/simulated comparison during objective evaluations (red points / blue curves, axes labelled “Chemical shift, Hz”). The top-level function has no declared output: it displays the final parameter vector, while the objective scalar stays local. No result spectrum, residual figure, or measured convergence status is supplied by the source.

## Distinction from the companion example

This is the maleate dataset and uses its own starting vector and the source's hours estimate; the common model structure does not make its guesses transferable fit results for fumarate.
