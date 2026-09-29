# examples/dnp_sol/cross_effect_freq_scan_3.m

- MATLAB implementation: [examples/dnp_sol/cross_effect_freq_scan_3.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/cross_effect_freq_scan_3.m)

## Purpose

A three-spin TOTAPOL cross-effect dynamic nuclear polarisation (DNP) frequency scan, set up to reproduce Figure 2b of the [Journal of Magnetic Resonance paper](https://doi.org/10.1016/j.jmr.2011.09.047). The source describes electron rotating-frame dynamics with Nottingham DNP relaxation theory, citing [the relaxation-theory paper](https://doi.org/10.1007/s00723-012-0367-0). It explicitly warns that calculated intensities differ because its relaxation model differs from the paper's and because the stated geometry and interaction amplitudes have minor inconsistencies.

## Model and setup

The no-argument function `cross_effect_freq_scan_3()` sets `sys.magnet=3.4` and uses two electron spins and one proton (`{'E','E','1H'}`). The electron Zeeman scalars are 2.0023193 and 1.9992887; the proton scalar is 0. The coordinates, in spin order, are `[0,0,0]`, `[12.80,0,0]`, and `[-3.12,0,3.12]`. The source does not state coordinate units. The basis is the full `sphten-liouv` basis (`approximation='none'`).

Relaxation is set to `{'nottingham'}`, with `rlx_keep='secular'`, `equilibrium='zero'`, and `temperature=10`. The literal relaxation values are `nott_r1e=1e2`, `nott_r2e=1e5`, `nott_r1n=0.1`, and `nott_r2n=1e3`; no units are annotated for these values in the script.

## Calculation and output

The sequence parameters select electron irradiation, set `mw_pwr=2*pi*100e3`, and construct `mw_frq=2*pi*linspace(-350,350,1e4)*1e6`. This variant samples -350 to 350 MHz at 10,000 offsets, compared with 50,000 in `cross_effect_freq_scan_2`; the second electron Zeeman scalar also differs. The source does not annotate a unit for `mw_pwr`. Detection is proton `Lz`; electron `Lx` and `Lz` operators are supplied for microwave and electron observables, respectively. It fixes `orientation=[0 0 0]`, uses `method='lvn-backs'`, requests `needs={'aniso_eq'}`, and sets `g_ref` to the first electron Zeeman scalar.

The calculation is `crystal(spin_system,@dnp_freq_scan,parameters,'esr')`. It stores the result in the local variable `answer` and plots `real(answer)` against microwave-frequency offset, with the proton `S_z` expectation value on the vertical axis. The function declares no output argument; the observable is presented in the figure. The source estimates calculation time as seconds. Running it requires Spinach and its `dnp_freq_scan` sequence callback.
