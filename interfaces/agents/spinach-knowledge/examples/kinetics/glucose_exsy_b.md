# examples/kinetics/glucose_exsy_b.m

- Signature: `glucose_exsy_b()`

## Purpose

A two-dimensional EXSY forward simulation for transmembrane exchange of 3,3-difluoroglucose. The source says the parameters come from a fitting example set, but this function does not fit them: it simulates the EXSY signal and compares it with a locally loaded experimental spectrum. The source estimates a calculation time of seconds.

## Spin, relaxation, and exchange model

Eight 19F spins form four two-spin chemical subsystems: alpha inside, alpha outside, beta inside, and beta outside. The source supplies shifts, scalar couplings, and coordinates for each subsystem. See the complete set in the [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/kinetics/glucose_exsy_b.m).

It sets `sys.magnet=9.3933`; no unit is written beside that assignment.

The four-subsystem chemical-rate matrix has alpha entries `[-0.3438 0.1550; 0.3438 -0.1550]` and beta entries `[-0.7995 0.4350; 0.7995 -0.4350]` in the corresponding diagonal blocks. Rate units are not annotated. The source passes `[3.8034; 0; 14.2442; 0]` to `equilibrate`; this is the starting vector, not a reported equilibrated measurement.

The basis is `sphten-liouv` with no approximation. The relaxation configuration combines Redfield and `t1_t2`, keeps secular terms, sets zero equilibrium, and assigns `tau_c` values 0.9601e-9 for inside and 0.5255e-9 for outside subsystems. The source does not annotate the `tau_c` units. It also sets all eight `r1_rates` entries to zero and each `r2_rates` entry to 34.0359; units are not stated. The mixing-time parameter is 0.5, with no unit annotation.

## EXSY acquisition and processing

The no-argument function prepares the system and basis, sets `rho0` to the chemical-subsystem 19F `Lz` state, and calls `liquid(spin_system,@noesy,parameters,'nmr')`. Sequence values are offset -46681, sweep `[8650 8650]`, `npoints=[512 1024]`, zero filling `[1024 1024]`, and `axis_units='ppm'`. The source does not annotate the offset or sweep units.

Cosine and sine signals receive squared-cosine apodisation in both dimensions. The code Fourier-transforms along F2, constructs `f1_cos - 1i*f1_sin`, transforms along F1, and plots the negative real spectrum. It loads `spec` from `glucose_expt_b.mat`, applies `atranspose`, then `keep_rank(...,25)` to prepare the experimental matrix.

## Output and limits

The figures show simulated and experimental spectra side by side and a pointwise difference histogram with displayed range [-300, 300]. `glucose_expt_b.mat` must be available on the MATLAB path/current working directory. This is not the fitting example mentioned in the source comment and it prints no fit score; this page therefore reports source inputs and code outputs, not a fitted result or an unrun simulation outcome.
