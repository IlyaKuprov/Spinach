# examples/fitting/fluoroalkanes/syn_difluorobutane.m

- MATLAB implementation: [examples/fitting/fluoroalkanes/syn_difluorobutane.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fitting/fluoroalkanes/syn_difluorobutane.m)

- Signature: `syn_difluorobutane()`

## Purpose and data convention

Fit the 1H spectrum of syn-2,3-difluorobutane across two proton intervals (CH and methyl, abbreviated ME in the input variable names). The model contains two 19F spins, but this example loads and fits proton data only. See [10.1021/acs.joc.4c00670](https://doi.org/10.1021/acs.joc.4c00670); the source estimates calculation time in hours.

The entry point loads `syn_dfb_proton.mat` with `ch_axis_hz`, `ch_expt_data`, `me_axis_hz`, and `me_expt_data`. It normalises their integrals to −2 and −6, then concatenates both spectra and frequency axes. The source labels those axes Hz.

## Model and fit

`fminsearch` minimises the squared norm of the real-spectrum residual over a 9-element vector: parameters 1–7 are grouped scalar couplings (H–F, H–H, and F–F), parameter 8 is the Gaussian apodisation argument, and parameter 9 scales the simulated spectrum. The initial guess is `[23.95 6.47 0.90 4.36 18.15 47.88 -11.61 13.63 1.7]`; optimiser settings include `MaxIter=5000`, unlimited function evaluations, a `DiffMinChange` of `1e-3`, and central finite differences.

The Spinach model uses eight `1H` and two `19F` spins at `sys.magnet=11.7464`, a Zeeman–Hilbert basis without approximation, and `S3` groups on proton indices 1–3 and 4–6. A single 1H liquid-state acquisition uses `offset=1500`, `sweep=1800`, 4096 points, 32768-point zero filling, and Hz axes. After apodisation and scaling, the simulated spectrum is reversed and interpolated onto the concatenated experimental axis with `pchip`.

## Entry point and visible result

Run `syn_difluorobutane()` with the MAT file available. The function has no declared return value; it displays the optimiser vector and plots experimental points against the simulated line, including a full-spectrum panel and a 645–700 Hz zoom. Source comments estimate hours of calculation. No fit outcome, parameter uncertainty, or acceptance threshold is stated in the source.
