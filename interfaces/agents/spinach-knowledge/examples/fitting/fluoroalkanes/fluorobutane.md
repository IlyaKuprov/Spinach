# examples/fitting/fluoroalkanes/fluorobutane.m

- Signature: `fluorobutane()`

## Purpose

Fit experimental ^1H and ^19F NMR spectra by optimising J-couplings and separate spectral scale factors. The source comments describe this as fitting the ^1H NMR spectrum of **2-fluoropentane**, although the function and data files are named `fluorobutane`. See https://doi.org/10.1021/acs.joc.4c00670 for further details. Calculation time: hours.

## Workflow

- Load the fluorine spectrum and two proton spectral intervals from `fluorobutane_fluorine.mat` and `fluorobutane_proton.mat`. Normalise each by its integrated intensity, with a factor of two for the second proton interval, then concatenate the proton intervals and axes.
- Starting from an 11-parameter guess, use `fminsearch` with `MaxIter=5000` and `MaxFunEvals=Inf` to minimise the summed squared real-spectrum residuals. Nine parameters specify J-couplings; `A` and `B` scale the simulated ^1H and ^19F spectra independently.
- In each objective evaluation, construct a Spinach system containing nine ^1H spins and one ^19F spin at 11.7464 T. The fixed chemical shifts are 0.982 ppm (far CH3), 1.333 ppm (near CH3), 4.6075 ppm (CHF), 1.5949 and 1.6890 ppm (CH2), and −173.184 ppm (^19F). Model the two methyl groups with `S3` permutation symmetry in the unapproximated `zeeman-hilb` formalism.
- Simulate separate ^1H and ^19F acquisitions with `liquid(...,@acquire,...,'nmr')`, Gaussian-apodise each FID with width 6.0, zero-fill and Fourier-transform it, generate a ppm axis, and interpolate the result onto the experimental axis using `pchip`. The ^1H acquisition uses a 2000 Hz sweep, 4096 points and 32768-point zero fill; the ^19F acquisition uses a 250 Hz sweep, 512 points and 2048-point zero fill.
- Set proton-spectrum indices `3770:3820` to zero in both experiment and simulation before calculating the residual. Plot experimental and simulated proton and fluorine regions during fitting, print the trial parameters, and display the optimised parameters on completion.