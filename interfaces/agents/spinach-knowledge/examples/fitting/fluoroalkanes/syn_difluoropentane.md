# examples/fitting/fluoroalkanes/syn_difluoropentane.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fitting/fluoroalkanes/syn_difluoropentane.m`
- Signature: `syn_difluoropentane()`
- Total lines: 193

## Purpose

Fitting of 1H NMR spectrum of syn-2,4-difluoropentane with respect to J-couplings. See our paper for further details: Calculation time: hours

## Physical / mathematical content

- Fitting examples. These files formulate parameter-estimation workflows in which simulated spectra or observables are matched to data, usually through nonlinear optimisation, residual construction, and physically constrained parameterisations.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The file also defines local helper function(s): `errfun()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Load experimental data; implemented by `load('syn_dfp_fluorine.mat','axis_ppm','spec'); axis_expt_f=axis_ppm; spec_expt_f=spec`.
- Lines 19-20: Normalise and shift the data; implemented by `spec_expt_f= 7*spec_expt_f/max(spec_expt_f)`.
- Lines 24-26: Set the guess; implemented by `guess=[ 6.2780 23.5067 5.1361 7.0537 24.9899 16.9264 48.0475 1.6017 -14.5413 0.7437 0.8771 0.6311]`.
- Lines 28-29: Set optimiser options; implemented by `options=optimset('Display','iter','MaxIter',5000,'MaxFunEvals',Inf)`.
- Lines 31-32: Get the figure going; implemented by `kfigure(); scale_figure([2.00 1.75])`.
- Lines 34-37: Run the optimisation; implemented by `answer=fminsearch(@(x)errfun(axis_expt_f, spec_expt_f, axis_expt_ha,spec_expt_ha, axis_expt_hb,spec_expt_hb,x),guess,options)`.
- Lines 39-40: Display the result; implemented by `disp(answer)`.

### Key state/data transformations

- Lines 15: computes `load('syn_dfp_fluorine.mat','axis_ppm','spec'); axis_expt_f` using `load('syn_dfp_fluorine.mat','axis_ppm','spec'); axis_expt_f=axis_ppm; spec_expt_f=spec`.
- Lines 16: computes `load('syn_dfp_proton_a.mat','axis_ppm','spec'); axis_expt_ha` using `load('syn_dfp_proton_a.mat','axis_ppm','spec'); axis_expt_ha=axis_ppm; spec_expt_ha=spec`.
- Lines 17: computes `load('syn_dfp_proton_b.mat','axis_ppm','spec'); axis_expt_hb` using `load('syn_dfp_proton_b.mat','axis_ppm','spec'); axis_expt_hb=axis_ppm; spec_expt_hb=spec`.
- Lines 20: computes `spec_expt_f` using `spec_expt_f= 7*spec_expt_f/max(spec_expt_f)`.
- Lines 21: computes `spec_expt_ha` using `spec_expt_ha= 4*spec_expt_ha/max(spec_expt_ha)-0.1`.
- Lines 22: computes `spec_expt_hb` using `spec_expt_hb=10*spec_expt_hb/max(spec_expt_hb)`.
- Lines 25-26: computes `guess` using `guess=[ 6.2780 23.5067 5.1361 7.0537 24.9899 16.9264 48.0475 1.6017 -14.5413 0.7437 0.8771 0.6311]`.
- Lines 29: computes `options` using `options=optimset('Display','iter','MaxIter',5000,'MaxFunEvals',Inf)`.
- Lines 35-37: computes `answer` using `answer=fminsearch(@(x)errfun(axis_expt_f, spec_expt_f, axis_expt_ha,spec_expt_ha, axis_expt_hb,spec_expt_hb,x),guess,options)`.

### Local helper functions

- Line 45: `errfun()` — `function err=errfun(axis_expt_f, spec_expt_f,`. Silence Spinach
  - Representative operation: `axis_expt_ha,spec_expt_ha, axis_expt_hb,spec_expt_hb,params)`.
  - Representative operation: `axis_expt_hb,spec_expt_hb,params)`.

## Implementation structure

- Fitting of 1H NMR spectrum of syn-2,4-difluoropentane with
- respect to J-couplings. See our paper for further details:
- Calculation time: hours
- Load experimental data
- Normalise and shift the data
- Set the guess
- Set optimiser options
- Get the figure going
- Run the optimisation
- Display the result
- Least squares error function
- Silence Spinach

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `optimset()`, `kfigure()`, `scale_figure()`, `fminsearch()`, `errfun()`, `params()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `spec_theo_f()`, `spec_theo_ha()`, `spec_theo_hb()`.
