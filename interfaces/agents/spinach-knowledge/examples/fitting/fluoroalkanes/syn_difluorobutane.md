# examples/fitting/fluoroalkanes/syn_difluorobutane.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fitting/fluoroalkanes/syn_difluorobutane.m`
- Signature: `syn_difluorobutane()`
- Total lines: 155

## Purpose

Fitting of 1H NMR spectrum of syn-2,3-difluorobutane with respect to J-couplings. See our paper for further details: Calculation time: hours

## Physical / mathematical content

- Fitting examples. These files formulate parameter-estimation workflows in which simulated spectra or observables are matched to data, usually through nonlinear optimisation, residual construction, and physically constrained parameterisations.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The file also defines local helper function(s): `errfun()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-16: Load experimental data; implemented by `load('syn_dfb_proton.mat','ch_axis_hz','ch_expt_data', 'me_axis_hz','me_expt_data')`.
- Lines 18-19: Normalize the data; implemented by `ch_expt_data=-2*ch_expt_data/trapz(ch_axis_hz,ch_expt_data)`.
- Lines 22-23: Concatentate spectral intervals; implemented by `expt_data=[ch_expt_data; me_expt_data]`.
- Lines 26-27: Set the guess; implemented by `guess=[23.95 6.47 0.90 4.36 18.15 47.88 -11.61 13.63 1.7]`.
- Lines 29-30: Get a figure going; implemented by `kfigure(); scale_figure([1.0 1.5])`.
- Lines 32-34: Set optimizer options; implemented by `options=optimset('Display','iter','MaxIter',5000,'MaxFunEvals',Inf, 'DiffMinChange',1e-3,'FinDiffType','central')`.
- Lines 36-37: Run the optimisation; implemented by `answer=fminsearch(@(x)errfun(axis_hz,expt_data,x),guess,options)`.
- Lines 39-40: Display the result; implemented by `disp(answer)`.

### Key state/data transformations

- Lines 19: computes `ch_expt_data` using `ch_expt_data=-2*ch_expt_data/trapz(ch_axis_hz,ch_expt_data)`.
- Lines 20: computes `me_expt_data` using `me_expt_data=-6*me_expt_data/trapz(me_axis_hz,me_expt_data)`.
- Lines 23: computes `expt_data` using `expt_data=[ch_expt_data; me_expt_data]`.
- Lines 24: computes `axis_hz` using `axis_hz=[ch_axis_hz; me_axis_hz]`.
- Lines 27: computes `guess` using `guess=[23.95 6.47 0.90 4.36 18.15 47.88 -11.61 13.63 1.7]`.
- Lines 33-34: computes `options` using `options=optimset('Display','iter','MaxIter',5000,'MaxFunEvals',Inf, 'DiffMinChange',1e-3,'FinDiffType','central')`.
- Lines 37: computes `answer` using `answer=fminsearch(@(x)errfun(axis_hz,expt_data,x),guess,options)`.

### Local helper functions

- Line 45: `errfun()` — `function err=errfun(axis_hz,expt_data,params)`. Silence Spinach
  - Representative operation: `sys.output='hush'`.
  - Representative operation: `sys.disable={'hygiene'}`.

## Implementation structure

- Fitting of 1H NMR spectrum of syn-2,3-difluorobutane with
- respect to J-couplings. See our paper for further details:
- Calculation time: hours
- Load experimental data
- Normalize the data
- Concatentate spectral intervals
- Set the guess
- Get a figure going
- Set optimizer options
- Run the optimisation
- Display the result
- Least squares error function

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `trapz()`, `kfigure()`, `scale_figure()`, `optimset()`, `fminsearch()`, `errfun()`, `params()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `sim_spec()`, `sweep2ticks()`.
