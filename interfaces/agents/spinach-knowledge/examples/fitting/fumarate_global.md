# examples/fitting/fumarate_global.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fitting/fumarate_global.m`
- Signature: `fumarate_global()`
- Total lines: 134

## Purpose

Simultaneous fitting of 1H and 13C NMR spectra of a slightly asymmetric fumarate diester. Calculation time: minutes

## Physical / mathematical content

- Fitting examples. These files formulate parameter-estimation workflows in which simulated spectra or observables are matched to data, usually through nonlinear optimisation, residual construction, and physically constrained parameterisations.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The file also defines local helper function(s): `errfun()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Load experimental data; implemented by `load('fumarate_1h.mat', 'axis_hz','data_expt'); axis_expt_h=axis_hz; spec_expt_h=data_expt`.
- Lines 16-17: Normalise the data; implemented by `spec_expt_h=spec_expt_h/max(spec_expt_h)`.
- Lines 20-22: Set the guess; implemented by `guess=[ -0.0325 0.0300 0.0040 -0.0030 70.9989 -2.7853 166.6826 15.6814 0.0366 0.0242 10.4768 5.4905]`.
- Lines 24-25: Set optimiser options; implemented by `options=optimset('Display','iter','MaxIter',5000,'MaxFunEvals',Inf)`.
- Lines 27-29: Run the optimisation; implemented by `answer=fminsearch(@(x)errfun(axis_expt_h,axis_expt_c, spec_expt_h,spec_expt_c,x),guess,options)`.
- Lines 31-32: Display the result; implemented by `disp(answer)`.

### Key state/data transformations

- Lines 13: computes `load('fumarate_1h.mat', 'axis_hz','data_expt'); axis_expt_h` using `load('fumarate_1h.mat', 'axis_hz','data_expt'); axis_expt_h=axis_hz; spec_expt_h=data_expt`.
- Lines 14: computes `load('fumarate_13c.mat','axis_hz','data_expt'); axis_expt_c` using `load('fumarate_13c.mat','axis_hz','data_expt'); axis_expt_c=axis_hz; spec_expt_c=data_expt`.
- Lines 17: computes `spec_expt_h` using `spec_expt_h=spec_expt_h/max(spec_expt_h)`.
- Lines 18: computes `spec_expt_c` using `spec_expt_c=spec_expt_c/max(spec_expt_c)`.
- Lines 21-22: computes `guess` using `guess=[ -0.0325 0.0300 0.0040 -0.0030 70.9989 -2.7853 166.6826 15.6814 0.0366 0.0242 10.4768 5.4905]`.
- Lines 25: computes `options` using `options=optimset('Display','iter','MaxIter',5000,'MaxFunEvals',Inf)`.
- Lines 28-29: computes `answer` using `answer=fminsearch(@(x)errfun(axis_expt_h,axis_expt_c, spec_expt_h,spec_expt_c,x),guess,options)`.

### Local helper functions

- Line 37: `errfun()` — `function err=errfun(axis_expt_h,axis_expt_c,`. Silence Spinach
  - Representative operation: `spec_expt_h,spec_expt_c,params)`.
  - Representative operation: `sys.output='hush'`.

## Implementation structure

- Simultaneous fitting of 1H and 13C NMR spectra of a slightly
- asymmetric fumarate diester.
- Calculation time: minutes
- Load experimental data
- Normalise the data
- Set the guess
- Set optimiser options
- Run the optimisation
- Display the result
- Least squares error function
- Silence Spinach
- Absorb parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `optimset()`, `fminsearch()`, `errfun()`, `params()`, `spin()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `sweep2ticks()`, `subplot()`, `kxlabel()`.
