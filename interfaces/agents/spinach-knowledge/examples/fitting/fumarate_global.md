# examples/fitting/fumarate_global.m

- MATLAB implementation: [examples/fitting/fumarate_global.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fitting/fumarate_global.m)

- Signature: `fumarate_global()`

## Purpose

Jointly fits the proton and carbon spectra of a slightly asymmetric fumarate diester. The source estimates minutes for the example; it does not report a fitted result.

## Data and fit vector

The function expects `fumarate_1h.mat` and `fumarate_13c.mat` in the MATLAB working directory. Each file supplies `axis_hz` and `data_expt`; the script keeps the axes and divides each spectrum by its own maximum before fitting. The single 12-element trial vector is shared by both channels:

| Entries | Role in the model |
| --- | --- |
| 1–2 | Carbon chemical shifts `c_shift_1`, `c_shift_2` |
| 3–4 | Proton chemical shifts `h_shift_1`, `h_shift_2` |
| 5–8 | `j_cc`, far and near `j_ch`, and `j_hh` scalar-coupling entries |
| 9–10 | Independent proton and carbon signal scales `a_h`, `a_c` |
| 11–12 | Proton and carbon exponential-apodisation settings `lw_h`, `lw_c` |

Those names and roles follow the source assignments; it does not state units for the fitted shift, coupling, linewidth, or scale values. The starting vector is `[-0.0325, 0.0300, 0.0040, -0.0030, 70.9989, -2.7853, 166.6826, 15.6814, 0.0366, 0.0242, 10.4768, 5.4905]`; these are initial guesses, not reported estimates.

## Spin and acquisition model

The objective builds one four-spin `1H/13C/13C/1H` system at the field expression `2*pi*500.101412e6/spin('1H')`, using a full Zeeman-Hilbert basis (`approximation='none'`). It separately simulates proton and carbon liquid-state acquisitions with `@acquire`: each has zero offset, no decoupling, 600 sweep, 1024 acquired points, 4096-point zero fill, and an axis labelled in Hz. For each channel, the simulated FID is exponentially apodised, divided by 10, Fourier transformed, and multiplied by its fitted scale. The generated axes come from `sweep2ticks`; each simulated spectrum is then shape-preserving-cubic interpolated with `interp1(...,'pchip')` onto its corresponding experimental axis before plotting and residual evaluation. No unit for the trial-vector entries is inferred from the displayed-axis setting.

## Objective and what the entry point exposes

`fumarate_global()` calls `fminsearch` from the vector above, with iterative display, a 5000-iteration cap, and no finite function-evaluation cap. Its local objective is the sum of the squared Euclidean residual norms for the normalised proton and carbon spectra. It draws a two-panel comparison during objective evaluation (experimental points in red; simulated curves in blue) with both axes labelled “Chemical shift, Hz”. The top-level function declares no output argument: it displays the optimiser's final vector; the scalar objective is internal. The file contains no saved-fit or fit-quality result, and reaching a displayed endpoint alone is not a source-supported claim of convergence.
