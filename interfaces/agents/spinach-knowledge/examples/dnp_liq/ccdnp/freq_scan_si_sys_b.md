# examples/dnp_liq/ccdnp/freq_scan_si_sys_b.m

- MATLAB implementation: [examples/dnp_liq/ccdnp/freq_scan_si_sys_b.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_liq/ccdnp/freq_scan_si_sys_b.m)

## What it computes

A steady-state liquid-state DNP map for one proton coupled dipolarly to two exchange-coupled electrons. It scans microwave frequency offset and static-field magnitude; it is a simulation example, not a reader for experimental data.

## Running and setup

Call the no-argument MATLAB function `freq_scan_si_sys_b()` with Spinach on the MATLAB path. The model, scan grids, and controls are assigned in the function; there are no user parameters or input files. It calls Spinach's `create`, `basis`, `state`, `operator`, `equilibrium`, and `liquid` routines, with `dnp_freq_scan` supplied as the ESR callback; that callback must also be available on the MATLAB path. The field loop is written as `parfor`, so parallel execution depends on MATLAB's parallel-computing setup.

## Model and fixed settings

- Isotopes are `{'1H','E','E'}`. The proton Zeeman principal values are `[0 10 20]` with zero Euler angles. Electron Zeeman values are `[1.9778 1.9776 1.9776]` and `[2.0068 2.0038 2.0038]`; their Euler triples are `[-0.180 0.017 0.194]` and `[0.632 0.783 1.086]`.
- The scalar electron-electron coupling is assigned as `5e6`. Coordinates are `[0 0 0]`, `[6.000 0.030 0.317]`, and `[-6.000 -0.038 0.535]`. The source gives no units for these values or coordinates.
- Basis: `sphten-liouv` with `none` approximation. Relaxation is Redfield, equilibrium mode `zero`, and relaxation terms are kept in the `secular` setting. The source sets temperature to `298`, `tau_c={100e-12}` (commented as TEMPOL in water), and relaxation-integration tolerance to `1e-10`; it does not state units for these literals.
- The driven electron is `E`; `mw_pwr=2*pi*2e6`, method `lvn-backs`, and `rho_eq` is requested. The reference g value is the mean of electron 1's three listed principal values. The power literal has no unit stated in the file.
- Microwave offsets use `2*pi*linspace(-15,15,512)*1e6`; the plotted axis expresses this as MHz relative to the electron-1 reference. The field grid is 128 points from 1 to 30 Tesla (the source explicitly labels the field in Tesla).

## Output and interpretation

At each field, `liquid(...,'esr')` supplies one response per microwave offset. The script divides each field column by the equilibrium proton-coil reference `coil'*rho_eq`, then plots `real(answer)` as a field-versus-offset image with a colorbar labelled steady-state proton DNP. A useful source-specific clarification is that rows correspond to the 512 offsets and columns to the 128 field values, matching the y- and x-axes passed to `imagesc`. The script creates a figure but does not save the grid or export the figure.

## Reference and caveat

The source cites https://doi.org/10.1016/j.jmr.2021.106940. no numerical DNP values are reported here.
