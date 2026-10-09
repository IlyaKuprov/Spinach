# examples/dnp_liq/ccdnp/tau_scan_si_sys_b.m

- MATLAB implementation: [examples/dnp_liq/ccdnp/tau_scan_si_sys_b.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_liq/ccdnp/tau_scan_si_sys_b.m)

## What it computes

This no-argument MATLAB function maps the steady-state proton response of a liquid-state DNP model against microwave-frequency offset and rotational correlation time. The spin system is a proton plus two exchange-coupled electrons, with anisotropic electron Zeeman tensors and dipolar interactions defined through coordinates. The source cites DOI: [10.1016/j.jmr.2021.106940](https://doi.org/10.1016/j.jmr.2021.106940). Its comment says calculation time is seconds; this was not timed here.

## Setup and parameter set

Use MATLAB with Spinach on the path and make the callback `dnp_freq_scan` available. The function takes no arguments and assigns the system internally; its `parfor` loop spans correlation times.

- Field `sys.magnet=14.1` T; isotopes `{'1H','E','E'}`. Proton Zeeman eigenvalues are `[0 10 20]` with Euler array `[0 0 0]`.
- Electron 1 eigenvalues `[1.9778 1.9776 1.9776]`, Euler array `[-0.180 0.017 0.194]`; electron 2 eigenvalues `[2.0068 2.0038 2.0038]`, Euler array `[0.632 0.783 1.086]`. Units/conventions for the eigenvalue and Euler arrays are not stated in this file.
- Scalar coupling array is initialised empty except `inter.coupling.scalar{2,3}=5e6` (unit not stated). Coordinates are `[0 0 0]`, `[6.000 0.030 0.317]`, and `[-6.000 -0.038 0.535]`; coordinate units are not given.
- Basis `sphten-liouv`, approximation `none`; Redfield relaxation; equilibrium `zero`; `rlx_keep='secular'`; temperature value `298` (no unit annotated); relaxation-integration tolerance `1e-10`.
- Sequence settings: `parameters.spins={'E'}`, `mw_pwr=2*pi*1e6` (unit not annotated), `method='lvn-backs'`, `needs={'rho_eq'}`; `g_ref` is the mean of electron 1's Zeeman eigenvalues. Hygiene checks are disabled and Spinach output is hush.

## Scan and plotted quantity

The microwave offset vector is `2*pi*linspace(-15,15,512)*1e6`, shown in MHz. The correlation-time vector is `linspace(50e-12,500e-12,128)` seconds, shown as 50–500 ps. For each correlation time, the function creates the Spinach system and basis, forms a proton `Lz` coil state and electron `Lx/2` and `Lz` operators, evaluates `liquid(spin_system,@dnp_freq_scan,localpar,'esr')`, and normalises by the proton detection expectation in the isotropic `equilibrium(spin_system)`. It plots the real part as a heat map; the x-axis is `τ_c` in ps and the y-axis is microwave offset from `g_iso^(1)` in MHz. No output data or figure is saved by the function.

**Source-specific clarification:** unlike parameter set A, this set has a strongly separated second-electron Zeeman tensor (`[2.0068 2.0038 2.0038]`), a different coordinate geometry and scalar-coupling input, and a symmetric `-15 to +15 MHz` scan.

**Caveats:** the source provides a runtime estimate and resulting DNP surface, which are estimates rather than guarantees. Values whose units or angular conventions are not documented above are kept as source values rather than assigned inferred units.
