# examples/fitting/fluoroalkanes/anti_difluoroheptane.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fitting/fluoroalkanes/anti_difluoroheptane.m`
- Signature: `anti_difluoroheptane()`
- Total lines: 208

## Purpose

Fitting of 1H NMR spectrum of anti-3,5-difluoroheptane with respect to J-couplings. See our paper for further details: Methyl groups are ghosted out because they do not influence the signals in question. Calculation time: hours

## Physical / mathematical content

- Fitting examples. These files formulate parameter-estimation workflows in which simulated spectra or observables are matched to data, usually through nonlinear optimisation, residual construction, and physically constrained parameterisations.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The file also defines local helper function(s): `errfun()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Load experimental data; implemented by `load('anti_dfh_fluorine.mat','axis_ppm','spec'); axis_expt_f =axis_ppm; spec_expt_f =spec`.
- Lines 22-23: Normalise the data; implemented by `spec_expt_f= spec_expt_f/max(spec_expt_f)`.
- Lines 27-29: Set the guess; implemented by `guess=[-14.4404 2.1802 10.0564 14.0478 38.0147 27.7945 18.4317 1.2295 3.1159 -15.1172 11.2721 49.6307 18.9158 4.5308 7.5739]`.
- Lines 31-32: Set optimiser options; implemented by `options=optimset('Display','iter','MaxIter',5000,'MaxFunEvals',Inf)`.
- Lines 34-35: Get the figure going; implemented by `kfigure(); scale_figure([2.00 1.75])`.
- Lines 37-40: Run the optimisation; implemented by `answer=fminsearch(@(x)errfun(axis_expt_f, spec_expt_f, axis_expt_ha,spec_expt_ha, axis_expt_hb,spec_expt_hb,x),guess,options)`.
- Lines 42-43: Display the result; implemented by `disp(answer)`.

### Key state/data transformations

- Lines 18: computes `load('anti_dfh_fluorine.mat','axis_ppm','spec'); axis_expt_f` using `load('anti_dfh_fluorine.mat','axis_ppm','spec'); axis_expt_f =axis_ppm; spec_expt_f =spec`.
- Lines 23: computes `spec_expt_f` using `spec_expt_f= spec_expt_f/max(spec_expt_f)`.
- Lines 24: computes `spec_expt_ha` using `spec_expt_ha=spec_expt_ha/max(spec_expt_ha)`.
- Lines 25: computes `spec_expt_hb` using `spec_expt_hb=spec_expt_hb/max(spec_expt_hb)`.
- Lines 28-29: computes `guess` using `guess=[-14.4404 2.1802 10.0564 14.0478 38.0147 27.7945 18.4317 1.2295 3.1159 -15.1172 11.2721 49.6307 18.9158 4.5308 7.5739]`.
- Lines 32: computes `options` using `options=optimset('Display','iter','MaxIter',5000,'MaxFunEvals',Inf)`.
- Lines 38-40: computes `answer` using `answer=fminsearch(@(x)errfun(axis_expt_f, spec_expt_f, axis_expt_ha,spec_expt_ha, axis_expt_hb,spec_expt_hb,x),guess,options)`.

### Local helper functions

- Line 48: `errfun()` — `function err=errfun(axis_expt_f, spec_expt_f,`. Silence the output
  - Representative operation: `axis_expt_ha,spec_expt_ha, axis_expt_hb,spec_expt_hb,params)`.
  - Representative operation: `axis_expt_hb,spec_expt_hb,params)`.

## Implementation structure

- Fitting of 1H NMR spectrum of anti-3,5-difluoroheptane with
- respect to J-couplings. See our paper for further details:
- Methyl groups are ghosted out because they do not influence
- the signals in question.
- Calculation time: hours
- Load experimental data
- Normalise the data
- Set the guess
- Set optimiser options
- Get the figure going
- Run the optimisation
- Display the result

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `optimset()`, `kfigure()`, `scale_figure()`, `fminsearch()`, `errfun()`, `params()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `spec_theo_f()`, `spec_theo_ha()`, `spec_theo_hb()`.
