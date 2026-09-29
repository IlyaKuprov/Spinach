# examples/dnp_sol/solid_effect_freq_scan_1.m

- Signature: `solid_effect_freq_scan_1()`
- Source: [`examples/dnp_sol/solid_effect_freq_scan_1.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/solid_effect_freq_scan_1.m)

## Purpose

Runs a laboratory-frame steady-state DNP calculation for a fixed crystal orientation and plots 1H and 15N longitudinal signals over two microwave-frequency windows. The source describes a single 15N-labelled urea molecule at a specified orientation and distance from one electron, identifies the calculation as single-crystal (not powder) and estimates minutes for runtime.

## Spin model and truncation

The configured isotope list is `{'E','15N','1H','1H','15N','1H','1H'}`: one electron, two 15N nuclei, and four 1H nuclei. The source sets `sys.magnet=3.4`. Coordinates in isotope order are:

`[0.00000000  0.00000000 10.14358975]`
`[-0.07640311 1.16112702 -0.61556225]`
`[0.08533754 1.99241453 -0.06489225]`
`[0.38824423 1.16155815 -1.51333625]`
`[0.07640311 -1.16112702 -0.61556225]`
`[-0.08533754 -1.99241453 -0.06489225]`
`[-0.38824423 -1.16155815 -1.51333625]`

The basis is `sphten-liouv` with `IK-0`, inter-level limit `inter_level=4`, and projections `[-2 -1 0 +1 +2]`. Relaxation uses the secular Weizmann model, zero equilibrium, and `inter.temperature=4.2`. The parameters are `weiz_r1e=1e2`, `weiz_r1n=0.1`, `weiz_r2e=1e5`, `weiz_r2n=1e3`; both `weiz_r1d` and `weiz_r2d` are `1e-3*ones(7,7)`. The source describes relaxation as accounting for T1, T2, and dipolar relaxation; it does not annotate units for the rates or temperature.

## Frequency calculation

The sequence acts on `E`, with `mw_pwr=2*pi*100e3`, electron Lx and Lz operators, fixed orientation `[pi/4 pi/5 pi/6]`, method `lvn-backs`, prerequisite `aniso_eq`, and `g_ref=spin_system.tols.freeg`. Detection combines the 1H and 15N Lz states into two observable columns.

The frequency vector concatenates two 100-point ranges: `linspace(144.0,145.5,100)` and `linspace(14.0,15.5,100)`, each multiplied by `2*pi*1e6` for `parameters.mw_frq`. The source labels the plotted horizontal coordinates in MHz. The calculation calls `crystal(spin_system,@dnp_freq_scan,parameters,'esr')`; it is therefore evaluated at the one configured orientation, rather than averaged over an orientation grid.

## Result and scope

The function makes a 2-by-2 figure: rows correspond to the 1H and 15N signal columns, and columns show the 144.0–145.5 MHz and 14.0–15.5 MHz windows. Each panel plots the real part of the corresponding 100-row segment of `answer`. The no-argument function has no declared MATLAB return value; its visible result is the figure. Its result depends on the explicit seven-spin model, the basis/relaxation truncations, the selected orientation, and the two frequency windows.
