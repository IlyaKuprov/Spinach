# examples/dnp_liq/ccdnp/states_vs_tauc_si_sys_b.m

- MATLAB implementation: [examples/dnp_liq/ccdnp/states_vs_tauc_si_sys_b.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_liq/ccdnp/states_vs_tauc_si_sys_b.m)

- Signature: `states_vs_tauc_si_sys_b()` (no arguments)
- Reference: [Journal of Magnetic Resonance 2021, 106940](https://doi.org/10.1016/j.jmr.2021.106940)
- The source comments that calculation takes seconds; this is not a benchmark guarantee.

## What this SI variant calculates

This file runs the supplementary-information B parameterisation of the same three-spin liquid-state DNP steady-state sweep. It tracks 12 selected observables versus rotational correlation time for a proton and two exchanging electrons. Choose it when the B-system tensors, geometry, and fixed microwave setting are required; it is not an alias for the main-text or SI A model.

## Running and model inputs

Run in MATLAB with Spinach and its example/plotting functions on the MATLAB path. The zero-argument function embeds the model and sweep; no data file is read. The loop is written as MATLAB `parfor`, with no pool setup in the file.

- Field: 14.1 T; spin order: `{'1H','E','E'}`. Nuclear Zeeman values are `[0 10 20]` with Euler angles `[0 0 0]`.
- Electron 1 uses `[1.9778 1.9776 1.9776]` at `[-0.180 0.017 0.194]`; electron 2 uses `[2.0068 2.0038 2.0038]` at `[0.632 0.783 1.086]`.
- The scalar entry for the electron pair is `5e6`. Coordinates are `[0 0 0]`, `[6.000 0.030 0.317]`, and `[-6.000 -0.038 0.535]`. The source gives no units for the coupling, coordinates, or Euler angles.
- Use the full `sphten-liouv` basis without approximation. Relaxation settings are Redfield, zero equilibrium polarisation, secular retention, temperature entry 298, and `sys.tols.rlx_integration=1e-10` (commented as needing to be this tight); the file does not state the temperature unit.
- DNP setup: `spins={'E'}`, `mw_pwr=2*pi*1e6`, `mw_frq=2*pi*3.2e6`, `method='lvn-backs'`, `needs={'rho_eq'}`, and `g_ref` equal to the mean of electron 1's Zeeman entries. Units for these microwave expressions are not given in the source.

The 64-point correlation-time grid is `50e-12` to `500e-12` (the plot axis is ps). Each `parfor` iteration constructs the spin system/basis, sets the corresponding correlation time, and evaluates `liquid(spin_system,@dnp_freq_scan,locpar,'esr')`. The same 12-element coil-state grouping is plotted: electron transverse combinations, mixed proton/electron terms, and longitudinal terms.

## Output and limits

The three panels make the 12 plotted states explicit: panel 1 contains `E1+ + 2 E1+E2z`, `E1+ - 2 E1+E2z`, `E2+ + 2 E1zE2+`, and `E2+ - 2 E1zE2+`; panel 2 contains `2 NzE1+`, `2 NzE2+`, `4 NzE1+E2z`, and `4 NzE1zE2+`; panel 3 contains `2 NzE1z`, `2 NzE2z`, `4 NzE1zE2z`, and `Nz`. Each curve is the absolute value of its corresponding result entry, plotted as steady-state amplitude in arbitrary units against correlation time in ps. No result array is returned and no data file is saved. A distinguishing source detail is the markedly different second-electron tensor (beginning at 2.0068) alongside the first-electron values near 1.9776; do not substitute SI A values when reproducing SI B. The source provides a calculation-time comment of seconds, not a measured timing on a specified machine. Numerical curves have not been reproduced here.