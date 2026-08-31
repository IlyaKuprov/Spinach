# examples/fitting/fluoroalkanes/syn_difluoroheptane.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fitting/fluoroalkanes/syn_difluoroheptane.m`
- Signature: `syn_difluoroheptane()`
- Total lines: 208

## Purpose

Fitting of 1H NMR spectrum of syn-3,5-difluoroheptane with respect to J-couplings. See our paper for further details: Methyl groups are ghosted out because they do not influence the signals in question. Calculation time: hours

## Physical / mathematical content

- Fitting examples. These files formulate parameter-estimation workflows in which simulated spectra or observables are matched to data, usually through nonlinear optimisation, residual construction, and physically constrained parameterisations.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The file also defines local helper function(s): `errfun()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Load experimental data; implemented by `load('syn_dfh_fluorine.mat','axis_ppm','spec'); axis_expt_f=axis_ppm; spec_expt_f=spec`.
- Lines 22-23: Normalize the data; implemented by `spec_expt_f= spec_expt_f/max(spec_expt_f)`.
- Lines 27-29: Set the guess; implemented by `guess=[ 1.0552 1.1935 0.9043 -13.8188 4.8881 7.0776 4.4430 7.7255 -14.7896 18.3289 25.8358 48.5693 17.2030 30.4926 1.8459]`.
- Lines 31-32: Set optimiser options; implemented by `options=optimset('Display','iter','MaxIter',5000,'MaxFunEvals',Inf)`.
- Lines 34-35: Get the figure going; implemented by `kfigure(); scale_figure([2.00 1.75])`.
- Lines 37-40: Run the optimisation; implemented by `answer=fminsearch(@(x)errfun(axis_expt_f, spec_expt_f, axis_expt_ha,spec_expt_ha, axis_expt_hb,spec_expt_hb,x),guess,options)`.
- Lines 42-43: Display the result; implemented by `disp(answer)`.

### Key state/data transformations

- Lines 18: computes `load('syn_dfh_fluorine.mat','axis_ppm','spec'); axis_expt_f` using `load('syn_dfh_fluorine.mat','axis_ppm','spec'); axis_expt_f=axis_ppm; spec_expt_f=spec`.
- Lines 19: computes `load('syn_dfh_proton_a.mat','axis_ppm','spec'); axis_expt_ha` using `load('syn_dfh_proton_a.mat','axis_ppm','spec'); axis_expt_ha=axis_ppm; spec_expt_ha=spec`.
- Lines 20: computes `load('syn_dfh_proton_b.mat','axis_ppm','spec'); axis_expt_hb` using `load('syn_dfh_proton_b.mat','axis_ppm','spec'); axis_expt_hb=axis_ppm; spec_expt_hb=spec`.
- Lines 23: computes `spec_expt_f` using `spec_expt_f= spec_expt_f/max(spec_expt_f)`.
- Lines 24: computes `spec_expt_ha` using `spec_expt_ha=spec_expt_ha/max(spec_expt_ha)`.
- Lines 25: computes `spec_expt_hb` using `spec_expt_hb=spec_expt_hb/max(spec_expt_hb)`.
- Lines 28-29: computes `guess` using `guess=[ 1.0552 1.1935 0.9043 -13.8188 4.8881 7.0776 4.4430 7.7255 -14.7896 18.3289 25.8358 48.5693 17.2030 30.4926 1.8459]`.
- Lines 32: computes `options` using `options=optimset('Display','iter','MaxIter',5000,'MaxFunEvals',Inf)`.
- Lines 38-40: computes `answer` using `answer=fminsearch(@(x)errfun(axis_expt_f, spec_expt_f, axis_expt_ha,spec_expt_ha, axis_expt_hb,spec_expt_hb,x),guess,options)`.

### Local helper functions

- Line 48: `errfun()` — `function err=errfun(axis_expt_f,spec_expt_f,`. Silence the output
  - Representative operation: `axis_expt_ha,spec_expt_ha, axis_expt_hb,spec_expt_hb,params)`.
  - Representative operation: `axis_expt_hb,spec_expt_hb,params)`.

## Implementation structure

- Fitting of 1H NMR spectrum of syn-3,5-difluoroheptane with
- respect to J-couplings. See our paper for further details:
- Methyl groups are ghosted out because they do not influence
- the signals in question.
- Calculation time: hours
- Load experimental data
- Normalize the data
- Set the guess
- Set optimiser options
- Get the figure going
- Run the optimisation
- Display the result

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `optimset()`, `kfigure()`, `scale_figure()`, `fminsearch()`, `errfun()`, `params()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `spec_theo_f()`, `spec_theo_ha()`, `spec_theo_hb()`.
