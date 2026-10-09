# examples/fitting/fluoroalkanes/anti_difluorobutane.m

- MATLAB implementation: [examples/fitting/fluoroalkanes/anti_difluorobutane.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fitting/fluoroalkanes/anti_difluorobutane.m)

- Signature: `anti_difluorobutane()`

## Purpose

Fit the 1H NMR spectrum of anti-2,3-difluorobutane by varying scalar couplings. The source cites [the associated paper](https://doi.org/10.1021/acs.joc.4c00670); its runtime comment says hours.

## Experimental data and spin model

The function loads `anti_dfb_proton.mat`, which must provide `ch_axis_hz`, `ch_expt_data`, `me_axis_hz`, and `me_expt_data`. It normalises the two signed spectra by their integrals using factors `-2` and `-6`, concatenates their data and Hz axes, then convolves the concatenated data with the normalised kernel `exp(-10*linspace(-1,1,100).^2)` to broaden the Z1-shim problem noted in the source.

The model has eight 1H and two 19F spins at a field of `11.7464`. The six methyl protons have shift `1.32375`, the other two protons `4.61040`, and both fluorines `0.00`. Two S3 groups cover proton sets `[1 2 3]` and `[4 5 6]`; the basis uses the Zeeman-Hilbert formalism without approximation.

## Fit parameterisation

The initial vector is `[24.08 6.49 1.44 3.59 15.73 47.76 -13.58 26.56 1.7]`. Its first seven entries are, in order, the near methyl–F coupling, methyl–methine H coupling, far methyl–F coupling, three-bond H–H, three-bond F–H, two-bond F–H, and three-bond F–F couplings. Entries 8 and 9 are the Gaussian linewidth and spectral amplitude. `fminsearch` minimises the squared norm of the difference between the real experimental and simulated spectra (`MaxIter=5000`, `MaxFunEvals=Inf`).

## Simulation and output

Each objective evaluation constructs one 1H liquid-acquisition simulation at offset `1500`, sweep `1800`, with `4096` acquired points and zero filling to `32768`; the source declares Hz axes and reverses the axis direction. It applies Gaussian apodisation using the fitted linewidth, Fourier transforms and scales by the fitted amplitude, then uses `pchip` interpolation onto the concatenated experimental axis. The figure compares the experiment and simulation in the two frequency windows `[2255 2355]` and `[635 690]` Hz; the final parameter vector is displayed.

The entry point has no output argument. The source defines a local unconstrained search and a residual objective, but does not provide a measured fit result or parameter uncertainties. The two input regions form one concatenated fit against a single simulated spectrum, rather than two independently optimised experiments.
