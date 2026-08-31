# examples/fitting/fluoroalkanes/fluorobutane.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fitting/fluoroalkanes/fluorobutane.m`
- Signature: `fluorobutane()`
- Total lines: 196

## Purpose

Fitting of 1H NMR spectrum of 2-fluoropentane with respect to J-couplings. See our paper for further details: Calculation time: hours

## Physical / mathematical content

- Fitting examples. These files formulate parameter-estimation workflows in which simulated spectra or observables are matched to data, usually through nonlinear optimisation, residual construction, and physically constrained parameterisations.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The file also defines local helper function(s): `errfun()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Load experimental data; implemented by `load('fluorobutane_fluorine.mat','spec_f','axis_f')`.
- Lines 19-20: Normalise the data; implemented by `spec_ch= -1*spec_ch/trapz(axis_ch,spec_ch)`.
- Lines 24-25: Concatenate spectral intervals; implemented by `spec_h=[spec_ch; spec_ch2]`.
- Lines 28-30: Set the guess; implemented by `guess=[23.9529 6.2219 7.4903 7.1692 4.9608 17.4939 26.2802 48.6869 -14.0924 1.8099 4.3571 ]`.
- Lines 32-33: Set optimiser options; implemented by `options=optimset('Display','iter','MaxIter',5000,'MaxFunEvals',Inf)`.
- Lines 35-36: Get the figure going; implemented by `kfigure(); scale_figure([2.00 1.75])`.
- Lines 38-39: Run the optimisation; implemented by `answer=fminsearch(@(x)errfun(axis_h,spec_h,axis_f,spec_f,x),guess,options)`.
- Lines 41-42: Display the result; implemented by `disp(answer)`.

### Key state/data transformations

- Lines 20: computes `spec_ch` using `spec_ch= -1*spec_ch/trapz(axis_ch,spec_ch)`.
- Lines 21: computes `spec_ch2` using `spec_ch2=-2*spec_ch2/trapz(axis_ch2,spec_ch2)`.
- Lines 22: computes `spec_f` using `spec_f= -1*spec_f/trapz(axis_f,spec_f)`.
- Lines 25: computes `spec_h` using `spec_h=[spec_ch; spec_ch2]`.
- Lines 26: computes `axis_h` using `axis_h=[axis_ch; axis_ch2]`.
- Lines 29-30: computes `guess` using `guess=[23.9529 6.2219 7.4903 7.1692 4.9608 17.4939 26.2802 48.6869 -14.0924 1.8099 4.3571 ]`.
- Lines 33: computes `options` using `options=optimset('Display','iter','MaxIter',5000,'MaxFunEvals',Inf)`.
- Lines 39: computes `answer` using `answer=fminsearch(@(x)errfun(axis_h,spec_h,axis_f,spec_f,x),guess,options)`.

### Local helper functions

- Line 47: `errfun()` — `function err=errfun(axis_h,expt_h,axis_f,expt_f,params)`. Silence Spinach
  - Representative operation: `sys.output='hush'`.
  - Representative operation: `sys.disable={'hygiene'}`.

## Implementation structure

- Fitting of 1H NMR spectrum of 2-fluoropentane with respect
- to J-couplings. See our paper for further details:
- Calculation time: hours
- Load experimental data
- Normalise the data
- Concatenate spectral intervals
- Set the guess
- Set optimiser options
- Get the figure going
- Run the optimisation
- Display the result
- Least squares error function

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `trapz()`, `optimset()`, `kfigure()`, `scale_figure()`, `fminsearch()`, `errfun()`, `params()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `spin()`, `sim_h()`.
