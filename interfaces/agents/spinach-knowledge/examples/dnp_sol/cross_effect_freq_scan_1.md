# examples/dnp_sol/cross_effect_freq_scan_1.m

- MATLAB implementation: [examples/dnp_sol/cross_effect_freq_scan_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/cross_effect_freq_scan_1.m)

- **Call:** `cross_effect_freq_scan_1()` (zero input arguments; the function plots the result and returns no explicit value).

## Purpose and source context

This TOTAPOL-based cross-effect DNP example calculates the proton response during electron rotating-frame irradiation using Nottingham DNP relaxation theory. It is configured to reproduce Fig. 2c of [the cited *Journal of Magnetic Resonance* paper](https://doi.org/10.1016/j.jmr.2011.09.047). The source cautions that intensity differences arise from a different relaxation model and minor inconsistencies between the stated geometry and interaction amplitudes in the original paper. The relaxation theory reference is [the Nottingham DNP paper](https://doi.org/10.1007/s00723-012-0367-0). The source estimates seconds of calculation time.

## Model and relaxation

The assigned system field is `sys.magnet=3.4` T. The isotopes are `{'E','E','1H'}`; scalar Zeeman entries are `{2.0023193,2.0021091,0.0000000}` (dimensionless electron g factors and a proton shift in ppm). Cartesian coordinates are `[0,0,0]`, `[12.80,0,0]` and `[-3.12,0,3.12]` in ångström. The full `sphten-liouv` basis uses no approximation.

Relaxation is `nottingham` with `rlx_keep='secular'` and `equilibrium='zero'`. The Nottingham rates `nott_r1e=1e2`, `nott_r2e=1e5`, `nott_r1n=0.1`, and `nott_r2n=1e3` are in hertz; `temperature=10` is in kelvin.

## Frequency scan and output

The irradiation spin is `'E'`, with `mw_pwr=2*pi*100e3` and `mw_frq=2*pi*linspace(-350,350,5e4)*1e6`. The scan therefore has 50,000 values over the plotted range `-350` to `+350 MHz`; the source's frequency assignment includes the `2*pi` factor. Detection is `coil_state(spin_system,'Lz','1H','exact')`; the microwave and electron-Zeeman operators are `operator(spin_system,'Lx','E')` and `operator(spin_system,'Lz','E')`. The fixed orientation is `[0 0 0]`; method is `lvn-backs`, `g_ref` is the first electron Zeeman scalar, and `needs={'aniso_eq'}`. The function calls `crystal(spin_system,@dnp_freq_scan,parameters,'esr')` and plots `real(answer)` against frequency offset.

**Observable-label note:** as in the source, the receiver is set to `Lz` on `1H`, whereas the plotted axis is labelled as an $S_z$ expectation on $^1$H; these are retained as distinct source details, not reconciled here.

**Dependencies:** Spinach system/basis/relaxation/state/operator and crystal routines; the `dnp_freq_scan` sequence and plotting helpers `kfigure`, `kgrid`, `kxlabel`, and `kylabel`.
