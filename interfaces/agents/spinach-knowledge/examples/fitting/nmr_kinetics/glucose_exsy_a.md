# examples/fitting/nmr_kinetics/glucose_exsy_a.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fitting/nmr_kinetics/glucose_exsy_a.m`
- Signature: `glucose_exsy_a()`
- Total lines: 201

## Purpose

Fitting of 2,2,3,3-tetrafluoroglucose NOESY with respect to the reaction rates in a chemical exchange problem and the rotational correlation time within Redfield theory. Calculation time: minutes (iteration count is limited in this example file)

## Physical / mathematical content

- Fitting examples. These files formulate parameter-estimation workflows in which simulated spectra or observables are matched to data, usually through nonlinear optimisation, residual construction, and physically constrained parameterisations.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The file also defines local helper function(s): `exsy_err()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Get a figure going; implemented by `kfigure(); scale_figure([1.75 1.00])`.
- Lines 17-23: Set the initial guess; implemented by `guess=[0.9448 1.6667 -8.4603 -9.0902 1.0074 0.1052 0.1032 0.0787 0.1125 0.1076 0.1273 0.0636 0.1367 0.0952 0.1171 0.0818 0.1226 0.1223 0.1079 0.1548 0.0500 -5.8355 -17.7…`.
- Lines 25-27: Set optimiser options; implemented by `options=optimset('Display','iter','MaxIter',10, 'MaxFunEvals',Inf,'UseParallel',true)`.
- Lines 29-30: Run the optimisation; implemented by `answer=fminsearch(@exsy_err,guess,options)`.
- Lines 32-33: Display the result; implemented by `disp(answer)`.

### Key state/data transformations

- Lines 18-23: computes `guess` using `guess=[0.9448 1.6667 -8.4603 -9.0902 1.0074 0.1052 0.1032 0.0787 0.1125 0.1076 0.1273 0.0636 0.1367 0.0952 0.1171 0.0818 0.1226 0.1223 0.1079 0.1548 0.0500 -5.8355 -17.7…`.
- Lines 26-27: computes `options` using `options=optimset('Display','iter','MaxIter',10, 'MaxFunEvals',Inf,'UseParallel',true)`.
- Lines 30: computes `answer` using `answer=fminsearch(@exsy_err,guess,options)`.

### Local helper functions

- Line 37: `exsy_err()` — `function err=exsy_err(params)`. Hush up Spinach
  - Representative operation: `sys.output='hush'`.
  - Representative operation: `sys.disable={'hygiene'}`.

## Implementation structure

- Fitting of 2,2,3,3-tetrafluoroglucose NOESY with respect to the
- reaction rates in a chemical exchange problem and the rotational
- correlation time within Redfield theory.
- Calculation time: minutes (iteration count is limited
- in this example file)
- Get a figure going
- Set the initial guess
- Set optimiser options
- Run the optimisation
- Display the result
- Hush up Spinach
- Magnet field

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `kfigure()`, `scale_figure()`, `optimset()`, `fminsearch()`, `exsy_err()`, `num2cell()`, `params()`, `equilibrate()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `load()`, `rot90()`.
