# examples/dnp_liq/ccdnp/tau_scan_si_sys_a.m

- MATLAB implementation: [examples/dnp_liq/ccdnp/tau_scan_si_sys_a.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_liq/ccdnp/tau_scan_si_sys_a.m)

## What it computes

This no-argument MATLAB function makes a two-axis liquid-state DNP scan for a proton coupled to two exchange-coupled electrons: the steady-state proton response is evaluated over microwave-frequency offset and rotational correlation time. The source cites DOI: [10.1016/j.jmr.2021.106940](https://doi.org/10.1016/j.jmr.2021.106940). Its comment estimates calculation time as seconds; that is source guidance, not a reproduced timing.

## Setup and spin model

Run with MATLAB and Spinach on the MATLAB path, from an environment where the example callback `dnp_freq_scan` is also resolvable. The scan uses `parfor` across correlation times. There are no function arguments or file inputs; all parameters are assigned in the function.

- Field: `sys.magnet=14.1` T; isotopes: `{'1H','E','E'}`.
- Proton Zeeman eigenvalues `[0 10 20]` and Euler array `[0 0 0]`. Electron 1 eigenvalues `[1.977873 1.977798 1.977792]`, Euler array `[0 0 0]`; electron 2 eigenvalues `[1.977919 1.978000 1.978000]`, Euler array `[-0.59 -0.10 0.49]`. Units/conventions for these arrays are not specified in this file.
- The only nonzero scalar-coupling entry is `inter.coupling.scalar{2,3}=6.2e6` (unit not stated). The three `inter.coordinates` rows are `[0 0 0]`, `[7.03 0.0187 0.9820]`, and `[-7.03 0.2051 -1.0001]`; coordinate units are not stated.
- Basis: `sphten-liouv`, approximation `none`. Relaxation is `redfield`, equilibrium setting `zero`, retained relaxation terms `secular`, and `inter.temperature=298` (the source gives no unit). The relaxation-integration tolerance is `1e-10`; the source comments that it needs to be this tight.
- Sequence settings: electron channel `{'E'}`, `parameters.mw_pwr=2*pi*1e6`, method `lvn-backs`, needs `{'rho_eq'}`, and reference g value is the mean of electron 1's Zeeman eigenvalues. The source does not annotate a unit for `mw_pwr`.

## Scan, normalisation, and output

The frequency-offset vector is `2*pi*linspace(-10,30,512)*1e6`; the figure expresses it in MHz. The correlation-time vector is `linspace(50e-12,500e-12,128)` seconds, plotted as 50–500 ps. For each time, the function builds the Spinach system and basis, detects proton `Lz`, defines electron `Lx/2` and `Lz` operators, calls `liquid(...,@dnp_freq_scan,...,'esr')`, then divides the response by the proton detection expectation in `equilibrium(spin_system)`. It stores the scan as a complex 512×128 `answer` array and plots `real(answer)` with `τ_c` (ps) horizontal and microwave offset relative to `g_iso^(1)` (MHz) vertical. The function creates a figure but does not save the numeric array or figure to a file. It disables the `hygiene` check and sets output to `hush`.

**Source-specific clarification:** this is the A parameter set, not a generic DNP model: it pairs two near-1.978 electron tensors and uses the asymmetric `-10 to +30 MHz` offset window, unlike the companion B scan.

**Caveats:** the example suppresses hygiene checks and uses hush output.  The listed Euler/coordinate/scalar-coupling conventions must be checked against the surrounding Spinach model if changing the parameter set.