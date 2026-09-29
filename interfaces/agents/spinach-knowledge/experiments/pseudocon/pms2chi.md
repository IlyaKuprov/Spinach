# experiments/pseudocon/pms2chi.m

Source: [MATLAB implementation](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/pseudocon/pms2chi.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=pms2chi.m)

## Purpose

Fits the magnetic-susceptibility tensor to observed paramagnetic shifts containing both contact and pseudocontact contributions, using Equation 10 of the cited model. Diamagnetic chemical shifts are not part of the fitted observations.

## Inputs and units

- `hfcs`: cell array of real symmetric 3-by-3 hyperfine tensors in Gauss, one per observation. The source notes that Gauss avoids dependence on the electron g-tensor and requires normalisation per unpaired electron in the S*A*I spin-Hamiltonian convention (as returned by `gparse`).
- `shifts`: real numeric vector of observed paramagnetic (contact plus pseudocontact) shifts in ppm.
- `isotopes`: cell array of character strings, one isotope label per observation (for example, `'13C'`).

The three arrays must have matching element counts. Each hyperfine tensor must be real, symmetric, and 3-by-3; isotope entries must be character strings.

## Fit and output

The objective is the sum of squared differences between each observed shift and `hfc2pms(hfcs{n},chi,isotopes{n})`. `fminunc` uses quasi-Newton/BFGS, starts six independent parameters at zero, allows at most 100 iterations and unlimited function evaluations, displays iterations, and enables parallel evaluation. The code imposes no physical bounds or regularisation.

The fitted `chi` is a general symmetric 3-by-3 tensor parameterised by six independent entries, so unlike `pcs2chi` it is not constrained to be traceless. It is returned in cubic Angstroms; `err` is the least-squares sum of squares.

## Reference

The source identifies Equation 10 in [DOI 10.1039/C4CP03106G](https://doi.org/10.1039/C4CP03106G).