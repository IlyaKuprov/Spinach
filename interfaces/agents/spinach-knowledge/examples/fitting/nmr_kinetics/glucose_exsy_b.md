# examples/fitting/nmr_kinetics/glucose_exsy_b.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fitting/nmr_kinetics/glucose_exsy_b.m`
- Signature: `glucose_exsy_b()`
- Total lines: 172

## Purpose

Fitting of 3,3-difluoroglucose NOESY with respect to the reaction rates in a chemical exchange and the rotational correlation times within Redfield theory. Calculation time: hours (iteration count is limited in this example file)

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
- Lines 17-21: Set the initial guess; implemented by `guess=[ 0.2331 0.1114 6.8446 0.1210 0.1423 0.1205 0.1697 0.0755 0.1158 0.0936 0.1256 26.9142 18.7751 25.5476 0.7948 0.4331 -9.0277 -9.2828 1.1924 32.9179]`.
- Lines 23-25: Set optimiser options; implemented by `options=optimset('Display','iter','MaxIter',10, 'MaxFunEvals',Inf,'UseParallel',true)`.
- Lines 27-28: Run the optimisation; implemented by `answer=fminsearch(@exsy_err,guess,options)`.
- Lines 30-31: Display the result; implemented by `disp(answer)`.
- Lines 33-34: Save figure; implemented by `savefig(gcf,'glucose_exsy_b.fig')`.

### Key state/data transformations

- Lines 18-21: computes `guess` using `guess=[ 0.2331 0.1114 6.8446 0.1210 0.1423 0.1205 0.1697 0.0755 0.1158 0.0936 0.1256 26.9142 18.7751 25.5476 0.7948 0.4331 -9.0277 -9.2828 1.1924 32.9179]`.
- Lines 24-25: computes `options` using `options=optimset('Display','iter','MaxIter',10, 'MaxFunEvals',Inf,'UseParallel',true)`.
- Lines 28: computes `answer` using `answer=fminsearch(@exsy_err,guess,options)`.

### Local helper functions

- Line 38: `exsy_err()` — `function err=exsy_err(params)`. Hush up Spinach
  - Representative operation: `sys.output='hush'`.
  - Representative operation: `sys.disable={'hygiene'}`.

## Implementation structure

- Fitting of 3,3-difluoroglucose NOESY with respect to the reaction
- rates in a chemical exchange and the rotational correlation times
- within Redfield theory.
- Calculation time: hours (iteration count is limited
- in this example file)
- Get a figure going
- Set the initial guess
- Set optimiser options
- Run the optimisation
- Display the result
- Save figure
- Hush up Spinach

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `kfigure()`, `scale_figure()`, `optimset()`, `fminsearch()`, `savefig()`, `exsy_err()`, `num2cell()`, `params()`, `equilibrate()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `load()`.
