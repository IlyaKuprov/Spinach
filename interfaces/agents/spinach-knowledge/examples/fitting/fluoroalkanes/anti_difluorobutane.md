# examples/fitting/fluoroalkanes/anti_difluorobutane.m

- Signature: `anti_difluorobutane()`

## Purpose

Fit the 1H NMR spectrum of anti-2,3-difluorobutane with respect to J-couplings. See the paper for further details: https://doi.org/10.1021/acs.joc.4c00670

Calculation time: hours.

## Physical / mathematical content

- The spin system contains eight 1H and two 19F spins at a magnetic field of 11.7464. The six methyl protons have chemical shift 1.32375, the two other protons 4.61040, and the fluorines 0.00. Two S3 symmetry groups cover the methyl proton sets.
- Nine fitted parameters specify near and far methyl–fluorine couplings, methyl–methine proton coupling, three-bond H–H and F–H couplings, two-bond F–H coupling, three-bond F–F coupling, Gaussian linewidth, and spectral amplitude.

## Numerical / algorithmic content

- Load the CH and methyl spectral intervals from `anti_dfb_proton.mat`; normalise them by their integrals with factors -2 and -6, concatenate them, and convolve the data with a Gaussian-shaped filter to broaden out the Z1 shim problem.
- Start from `[24.08  6.49  1.44  3.59  15.73  47.76  -13.58  26.56  1.7]`. `fminsearch` minimises the squared norm of the difference between the real experimental and simulated spectra, with `MaxIter` set to 5000 and `MaxFunEvals` set to `Inf`.
- Simulate 1H acquisition using an offset of 1500 Hz, sweep of 1800 Hz, 4096 points, and zero filling to 32768 points. Apply Gaussian apodisation using the fitted linewidth, Fourier-transform and scale by the fitted amplitude, then interpolate onto the experimental frequency axis with `pchip`.
- Plot experiment and simulation in the frequency windows [2255 2355] and [635 690] Hz, and display the fitted parameters.