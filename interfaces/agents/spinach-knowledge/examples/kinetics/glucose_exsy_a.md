# examples/kinetics/glucose_exsy_a.m

- Signature: `glucose_exsy_a()`

## Purpose

A two-dimensional EXSY forward simulation for transmembrane exchange of 2,2,3,3-tetrafluoroglucose. The source says the parameters come from a fitting example set, but this function does not perform that fit: it constructs the spin/exchange model, simulates an EXSY signal, and compares it with a locally loaded experimental spectrum. The source estimates a calculation time of seconds.

## Spin, relaxation, and exchange model

Sixteen 19F spins form four four-spin chemical subsystems: alpha inside, alpha outside, beta inside, and beta outside. The source supplies per-spin chemical shifts, scalar couplings, and coordinates for those subsystems. See the full numeric definitions in the [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/kinetics/glucose_exsy_a.m).

It sets `sys.magnet=9.4`; no unit is written beside that assignment.

The four-subsystem chemical-rate matrix has alpha entries `[-1.0045 1.7738; 1.0045 -1.7738]` and beta entries `[-0.9304 1.4586; 0.9304 -1.4586]` in the corresponding diagonal blocks. The source does not annotate rate units. It passes `[3.2258; 0; 3.1902; 0]` to `equilibrate` as the starting concentration vector; the resulting concentration vector is computed at runtime and is not stated here as an observed value.

The basis is `sphten-liouv` with no approximation. Relaxation combines Redfield theory with secular retention, zero equilibrium, and `tau_c` values 4.137e-9 for the inside subsystems and 0.951e-9 for the outside subsystems (the source does not annotate the `tau_c` units). The mixing-time parameter is 0.5; its unit is not stated in the source.

## EXSY acquisition and processing

The callable function has no arguments. It creates the system and basis, sets `rho0` to the chemical-subsystem 19F `Lz` state, and calls `liquid(spin_system,@noesy,parameters,'nmr')`. It uses offset -49000, sweep `[8000 8000]`, `npoints=[256 256]`, zero filling `[1024 512]`, and `axis_units='ppm'`. The source does not state units for the offset or sweep assignments.

Both cosine and sine signal components receive squared-cosine apodisation in both dimensions. The code Fourier-transforms the components along F2, forms the States signal as `f1_cos - 1i*f1_sin`, transforms along F1, and plots the negative real spectrum. It loads `Expression1` from `glucose_expt_a.mat`, rotates the matrix by 180 degrees, divides it by 5, and applies `keep_rank(...,50)` before comparison.

## Output and limits

One figure places the simulated and experimental 2D spectra side by side. A second figure plots the pointwise difference histogram, with the displayed horizontal range fixed to [-300, 300]. This script needs `glucose_expt_a.mat` available on the MATLAB path/current working directory. It contains no fit procedure or printed fit score; no simulation outcome is quoted here.
