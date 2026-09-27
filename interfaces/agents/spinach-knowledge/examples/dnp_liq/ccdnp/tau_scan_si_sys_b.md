# examples/dnp_liq/ccdnp/tau_scan_si_sys_b.m

- Signature: `tau_scan_si_sys_b()`
- Reference: [Journal of Magnetic Resonance 2021, 106940](https://doi.org/10.1016/j.jmr.2021.106940)
- Calculation time: seconds

## Purpose

Computes and plots steady-state (^1mathrm{H}) DNP versus microwave-frequency offset and rotational correlation time for a three-spin system (one proton and two electrons). This system-B variant specifies its own electron Zeeman tensors, exchange coupling, and coordinates.

## Physical model

- The two electrons have anisotropic Zeeman tensors and are scalar exchange-coupled at (5	imes10^6); the listed coordinates define their geometry relative to the proton for the anisotropic hyperfine interactions.
- Relaxation is Redfield with secular retention, zero equilibrium polarization, and temperature 298 K; the relaxation integration tolerance is (10^{-10}).

## Calculation

The script sweeps 128 correlation times from 50 to 500 ps. At each value it builds the Spinach system and basis and runs a steady-state liquid-state DNP frequency scan using `dnp_freq_scan`. The 512 microwave offsets span -15 to 15 MHz. The proton response is normalized by its isotropic thermal-equilibrium reference, and the real part is plotted versus frequency offset and correlation time. The correlation-time sweep uses `parfor`.
