# examples/dnp_liq/ccdnp/states_vs_tauc_si_sys_a.m

- MATLAB implementation: [examples/dnp_liq/ccdnp/states_vs_tauc_si_sys_a.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_liq/ccdnp/states_vs_tauc_si_sys_a.m)

- Signature: `states_vs_tauc_si_sys_a()` (no arguments)
- Reference: [Journal of Magnetic Resonance 2021, 106940](https://doi.org/10.1016/j.jmr.2021.106940)
- The source comments that calculation takes seconds; this is not a benchmark guarantee.

## What this SI variant calculates

This is the supplementary-information A spin-system instance of the three-spin liquid-state DNP correlation-time sweep: one proton and two exchanging electrons, with electron-nuclear dipolar couplings, Redfield relaxation, and the same 12 selected steady-state spin observables as the main-text curve example. It is useful when the SI A parameterisation—not the main-text or SI B tensors—is the intended model.

## Running and model inputs

Run in MATLAB with Spinach and its example/plotting functions on the MATLAB path. The zero-argument function defines its model and sweep internally and reads no external data. It uses `parfor`; this file does not configure a parallel pool.

- Field is 14.1 T and the isotope order is `{'1H','E','E'}`. The nuclear Zeeman entries are `[0 10 20]` with zero Euler angles.
- Electron 1 has Zeeman values `[1.977873 1.977798 1.977792]` and zero Euler angles; electron 2 has `[1.977919 1.978000 1.978000]` and Euler angles `[-0.59 -0.10 0.49]`.
- The electron–electron scalar entry is `inter.coupling.scalar{2,3}=6.2e6`. Coordinates supplied for anisotropic hyperfine interactions are `[0 0 0]`, `[7.03 0.0187 0.9820]`, and `[-7.03 0.2051 -1.0001]`. The source does not state units for these values or the Euler angles.
- The basis is full `sphten-liouv` (`approximation='none'`). Relaxation is Redfield, with zero equilibrium polarisation, secular retention, temperature entry 298, and the explicitly tight integration tolerance `1e-10` (source comment: “Needs to be this tight”). No unit is attached to the temperature entry in the source.
- DNP parameters are `spins={'E'}`, `mw_pwr=2*pi*1e6`, `mw_frq=2*pi*15.4e6`, `method='lvn-backs'`, and `needs={'rho_eq'}`; `g_ref` is the mean of electron 1's Zeeman values. The microwave expressions have no unit stated in this file.

The 64 correlation times span `50e-12` to `500e-12` (displayed in ps). For each value, `parfor` rebuilds the spin system and basis and calls `liquid(spin_system,@dnp_freq_scan,locpar,'esr')`. The coil vector contains 12 observables: four electron-transverse terms, four nuclear/electron mixed terms, and four longitudinal terms.

## Output and limits

The three panels make the 12 plotted states explicit: panel 1 contains `E1+ + 2 E1+E2z`, `E1+ - 2 E1+E2z`, `E2+ + 2 E1zE2+`, and `E2+ - 2 E1zE2+`; panel 2 contains `2 NzE1+`, `2 NzE2+`, `4 NzE1+E2z`, and `4 NzE1zE2+`; panel 3 contains `2 NzE1z`, `2 NzE2z`, `4 NzE1zE2z`, and `Nz`. Each curve is the absolute value of its corresponding result entry, plotted as steady-state amplitude in arbitrary units against correlation time in ps. It neither returns the array nor writes numerical data to disk. This variant differs materially from the main-text model in its electron Zeeman values/orientations, exchange entry, coordinates, and fixed microwave expression. Only these configured points and observables are represented; the paper DOI is the source for broader context, not a claim that this page reproduces its results. Numerical curves have not been reproduced here.
