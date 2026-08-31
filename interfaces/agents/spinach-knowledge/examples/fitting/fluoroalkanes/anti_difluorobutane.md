# examples/fitting/fluoroalkanes/anti_difluorobutane.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fitting/fluoroalkanes/anti_difluorobutane.m`
- Signature: `anti_difluorobutane()`
- Total lines: 158

## Purpose

Fitting of 1H NMR spectrum of anti-2,3-difluorobutane with respect to J-couplings. See our paper for further details: Calculation time: hours

## Physical / mathematical content

- Fitting examples. These files formulate parameter-estimation workflows in which simulated spectra or observables are matched to data, usually through nonlinear optimisation, residual construction, and physically constrained parameterisations.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The file also defines local helper function(s): `errfun()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-16: Load experimental data; implemented by `load('anti_dfb_proton.mat','ch_axis_hz','ch_expt_data', 'me_axis_hz','me_expt_data')`.
- Lines 18-19: Normalise the data; implemented by `ch_expt_data=-2*ch_expt_data/trapz(ch_axis_hz,ch_expt_data)`.
- Lines 22-23: Concatenate spectral intervals; implemented by `expt_data=[ch_expt_data; me_expt_data]`.
- Lines 26-27: Broaden out Z1 shim problem; implemented by `filter=exp(-10*linspace(-1,1,100).^2)`.
- Lines 30-31: Set the guess; implemented by `guess=[24.08 6.49 1.44 3.59 15.73 47.76 -13.58 26.56 1.7]`.
- Lines 33-34: Get a figure going; implemented by `kfigure(); scale_figure([1.0 1.5])`.
- Lines 36-37: Set optimiser options; implemented by `options=optimset('Display','iter','MaxIter',5000,'MaxFunEvals',Inf)`.
- Lines 39-40: Run the optimisation; implemented by `answer=fminsearch(@(x)errfun(axis_hz,expt_data,x),guess,options)`.
- Lines 42-43: Display the result; implemented by `disp(answer)`.

### Key state/data transformations

- Lines 19: computes `ch_expt_data` using `ch_expt_data=-2*ch_expt_data/trapz(ch_axis_hz,ch_expt_data)`.
- Lines 20: computes `me_expt_data` using `me_expt_data=-6*me_expt_data/trapz(me_axis_hz,me_expt_data)`.
- Lines 23: computes `expt_data` using `expt_data=[ch_expt_data; me_expt_data]`.
- Lines 24: computes `axis_hz` using `axis_hz=[ch_axis_hz; me_axis_hz]`.
- Lines 27: computes `filter` using `filter=exp(-10*linspace(-1,1,100).^2)`.
- Lines 31: computes `guess` using `guess=[24.08 6.49 1.44 3.59 15.73 47.76 -13.58 26.56 1.7]`.
- Lines 37: computes `options` using `options=optimset('Display','iter','MaxIter',5000,'MaxFunEvals',Inf)`.
- Lines 40: computes `answer` using `answer=fminsearch(@(x)errfun(axis_hz,expt_data,x),guess,options)`.

### Local helper functions

- Line 48: `errfun()` — `function err=errfun(axis_hz,expt_data,params)`. Silence Spinach
  - Representative operation: `sys.output='hush'`.
  - Representative operation: `sys.disable={'hygiene'}`.

## Implementation structure

- Fitting of 1H NMR spectrum of anti-2,3-difluorobutane with
- respect to J-couplings. See our paper for further details:
- Calculation time: hours
- Load experimental data
- Normalise the data
- Concatenate spectral intervals
- Broaden out Z1 shim problem
- Set the guess
- Get a figure going
- Set optimiser options
- Run the optimisation
- Display the result

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `trapz()`, `kfigure()`, `scale_figure()`, `optimset()`, `fminsearch()`, `errfun()`, `params()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `sim_spec()`, `sweep2ticks()`.
