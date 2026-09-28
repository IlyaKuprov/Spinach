# examples/dnp_liq/ccdnp/tau_scan_main_text.m

- Signature: `tau_scan_main_text()`
- Reference: [Journal of Magnetic Resonance 2021, 106940](https://doi.org/10.1016/j.jmr.2021.106940)
- Calculation time: seconds (per source comment)

## Purpose

Computes steady-state proton DNP as a function of microwave-frequency offset and rotational correlation time for a three-spin system containing one proton and two electrons.

## Model and calculation

The example uses a 14.1 T field, anisotropic electron Zeeman tensors, scalar electron–electron exchange coupling of (3 × 10^6), and coordinates for the spin system. It uses a full sphten-liouv basis and Redfield relaxation with secular retention, zero equilibrium polarization, and temperature 298 K; the relaxation integration tolerance is (10^{-10}).

For each of 128 correlation times from 50 to 500 ps, a `parfor` loop builds the Spinach system and basis and runs a steady-state liquid-state frequency scan with `dnp_freq_scan`. The 512 microwave offsets span -5 to 10 MHz. The proton response is normalized by its isotropic thermal-equilibrium reference, and the real part is plotted against offset and correlation time.
