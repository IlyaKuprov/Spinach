# examples/fitting/fluoroalkanes/difluoropropane.m

- Signature: `difluoropropane()`

Fits the ^1H and ^19F NMR spectra of 1,3-difluoropropane with respect to J-couplings. See https://doi.org/10.1021/acs.joc.4c00670 for further details. Calculation time: hours.

The function loads a ^19F spectrum and two ^1H spectra from `difluoropropane_fluorine.mat`, `difluoropropane_proton_a.mat`, and `difluoropropane_proton_b.mat`, then normalises each to its maximum. It starts from the six-parameter guess `[0.3023 1.1605 0.7390 47.0061 5.7870 25.7705]` and runs `fminunc` with `MaxIter=5000` and `MaxFunEvals=Inf`. The first three parameters scale the simulated spectra; the last three set groups of J-couplings.

For each trial parameter set, `errfun` builds an eight-spin ^1H/^19F system at a magnetic field of 11.7464, using fixed chemical shifts and a `zeeman-hilb` basis with no approximation and two S2 spin symmetries. It simulates separate ^19F, ^1H A, and ^1H B acquisitions, applies exponential apodisation (8.0, 9.0, and 6, respectively), Fourier-transforms and reverses the spectra, and interpolates them onto the experimental ppm axes. The objective is the sum of squared spectral residuals across all three datasets. Each evaluation plots experimental points against simulated curves and displays the trial parameters; the function displays the optimised parameters when fitting finishes.