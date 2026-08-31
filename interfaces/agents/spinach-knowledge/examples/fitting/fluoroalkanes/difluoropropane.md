# examples/fitting/fluoroalkanes/difluoropropane.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fitting/fluoroalkanes/difluoropropane.m`
- Signature: `difluoropropane()`
- Total lines: 175

## Purpose

Fitting of 1H and 19F NMR spectrums of 1,3-difluoropropane with respect to J-couplings. See our paper for further details: Calculation time: hours

## Physical / mathematical content

- Fitting examples. These files formulate parameter-estimation workflows in which simulated spectra or observables are matched to data, usually through nonlinear optimisation, residual construction, and physically constrained parameterisations.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The file also defines local helper function(s): `errfun()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Load experimental data; implemented by `load('difluoropropane_fluorine.mat','axis_ppm','spec'); axis_expt_f=axis_ppm; spec_expt_f=spec`.
- Lines 19-20: Normalise the data; implemented by `spec_expt_f =spec_expt_f /max(spec_expt_f)`.
- Lines 24-25: Set the guess; implemented by `guess=[0.3023 1.1605 0.7390 47.0061 5.7870 25.7705]`.
- Lines 27-28: Set optimiser options; implemented by `options=optimset('Display','iter','MaxIter',5000,'MaxFunEvals',Inf)`.
- Lines 30-31: Get the figure going; implemented by `kfigure(); scale_figure([2.00 1.75])`.
- Lines 33-36: Run the optimisation; implemented by `answer=fminunc(@(x)errfun(axis_expt_f, spec_expt_f, axis_expt_ha,spec_expt_ha, axis_expt_hb,spec_expt_hb,x),guess,options)`.
- Lines 38-39: Display the result; implemented by `disp(answer)`.

### Key state/data transformations

- Lines 20: computes `spec_expt_f` using `spec_expt_f =spec_expt_f /max(spec_expt_f)`.
- Lines 21: computes `spec_expt_ha` using `spec_expt_ha=spec_expt_ha/max(spec_expt_ha)`.
- Lines 22: computes `spec_expt_hb` using `spec_expt_hb=spec_expt_hb/max(spec_expt_hb)`.
- Lines 25: computes `guess` using `guess=[0.3023 1.1605 0.7390 47.0061 5.7870 25.7705]`.
- Lines 28: computes `options` using `options=optimset('Display','iter','MaxIter',5000,'MaxFunEvals',Inf)`.
- Lines 34-36: computes `answer` using `answer=fminunc(@(x)errfun(axis_expt_f, spec_expt_f, axis_expt_ha,spec_expt_ha, axis_expt_hb,spec_expt_hb,x),guess,options)`.

### Local helper functions

- Line 44: `errfun()` — `function err=errfun(axis_expt_f, spec_expt_f,`. Silence the output
  - Representative operation: `axis_expt_ha,spec_expt_ha, axis_expt_hb,spec_expt_hb,params)`.
  - Representative operation: `axis_expt_hb,spec_expt_hb,params)`.

## Implementation structure

- Fitting of 1H and 19F NMR spectrums of 1,3-difluoropropane with
- respect to J-couplings. See our paper for further details:
- Calculation time: hours
- Load experimental data
- Normalise the data
- Set the guess
- Set optimiser options
- Get the figure going
- Run the optimisation
- Display the result
- Least squares error function
- Silence the output

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `optimset()`, `kfigure()`, `scale_figure()`, `fminunc()`, `errfun()`, `params()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `spec_theo_f()`, `spec_theo_ha()`, `spec_theo_hb()`.
