# experiments/pseudocon/pcs2chi.m

Source: [MATLAB implementation](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/pseudocon/pcs2chi.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=pcs2chi.m)

## Purpose

Fits the rank-2, traceless magnetic-susceptibility tensor to observed pseudocontact shifts using Equation 10 of the cited model. It is a least-squares fitting routine, not a pulse-sequence routine.

## Inputs and units

- `hfcs`: cell array of real symmetric 3-by-3 hyperfine tensors in Gauss, one per observation. The source notes that Gauss units avoid dependence on the electron g-tensor and that tensors should be normalised per unpaired electron, in the S*A*I spin-Hamiltonian convention (as returned by `gparse`).
- `shifts`: real numeric vector of pseudocontact shifts in ppm, excluding the diamagnetic contribution.
- `isotopes`: cell array of character strings naming the isotope for each observation, for example {'1H','13C'}.

The three inputs must contain the same number of observations. `hfcs` entries must be real symmetric 3-by-3 matrices; isotope entries must be character strings.

## Fit and output

For each tensor candidate, `hfc2pcs(hfcs{n},chi,isotopes{n})` supplies the predicted shift. The objective is the sum over observations of (observed shift - predicted shift)^2. `fminunc` uses the quasi-Newton algorithm with BFGS Hessian updates, starts all five independent parameters at zero, allows at most 100 iterations and unlimited function evaluations, and displays iterations. No physical bounds or regularisation are imposed.

The fitted `chi` is symmetric and traceless, parameterised as [x1 x2 x3; x2 x4 x5; x3 x5 -(x1+x4)]. Thus the routine returns only the anisotropic rank-2 part, in cubic Angstroms; the second output `err` is the least-squares sum of squares.

## Reference

The source identifies Equation 10 in [DOI 10.1039/C4CP03106G](https://doi.org/10.1039/C4CP03106G).