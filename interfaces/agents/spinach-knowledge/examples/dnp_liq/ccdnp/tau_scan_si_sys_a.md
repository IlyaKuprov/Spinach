# examples/dnp_liq/ccdnp/tau_scan_si_sys_a.m

- Signature: `tau_scan_si_sys_a()`
- Reference: [Concilio et al., Journal of Magnetic Resonance 2021, 106940](https://doi.org/10.1016/j.jmr.2021.106940)
- Calculation time: seconds

## Purpose

Computes and plots steady-state (^1mathrm{H}) DNP as a function of microwave-frequency offset and rotational correlation time for a three-spin system (one proton and two electrons). This variant uses the system-A electron Zeeman tensors and electron–electron exchange parameters from the cited example.

## Physical model

- The two electrons are exchange-coupled (scalar coupling set to (6.2 × 10^6)); the electron Zeeman tensors are anisotropic, and coordinates specify the electron–proton geometry used for the anisotropic hyperfine interactions.
- Relaxation uses Redfield theory with secular retention, zero equilibrium polarization, and temperature 298 K. The relaxation integration tolerance is set to (10^{-10}).

## Calculation

For each of 128 correlation times from 50 to 500 ps, the script constructs the Spinach system and basis, then runs a steady-state liquid-state DNP frequency scan with `dnp_freq_scan`. It uses 512 microwave offsets from -10 to 30 MHz, normalizes the proton signal by its isotropic thermal-equilibrium reference, and plots the real signal against offset and correlation time. The correlation-time sweep is parallelized with `parfor`.
