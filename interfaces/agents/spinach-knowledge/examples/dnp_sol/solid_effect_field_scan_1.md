# examples/dnp_sol/solid_effect_field_scan_1.m

- Signature: `solid_effect_field_scan_1()`
- Source: [`examples/dnp_sol/solid_effect_field_scan_1.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/solid_effect_field_scan_1.m)

## Purpose

Calculates a steady-state solid-effect DNP signal for a gadolinium electron–15N system over a magnetic-field-offset scan. The source describes the ordinate as the 15N longitudinal expectation value and estimates minutes for the run.

## Spin model and relaxation

The source sets `sys.magnet=9.4509` and isotopes `{'E8','15N'}`. The electron Zeeman principal values are all `1.9918`. It specifies an electron zero-field-splitting tensor with principal values `570e6*[-1/3 -1/3 2/3]` and zero Euler angles. The two coordinate rows are `[0 0 0]` and `[3 0 0]`; the source explicitly labels these coordinates in Angstrom.

The Liouville-space basis is `sphten-liouv`, with no approximation and the projection set `[-2 -1 0 1 2]`. Relaxation is configured as `t1_t2`, with `r1_rates={1e4 1e1}`, `r2_rates={1e7 1e3}`, `rlx_keep='diagonal'`, zero equilibrium, and `inter.temperature=40.2`. The source does not attach units to the temperature or rate values.

## Field-scan calculation

The sequence acts on `E8`; microwave parameters are `mw_pwr=1e5` and `mw_frq=-14e8` as assigned in the source. It detects the 15N Lz state and defines electron Lx and Lz operators as the microwave and electron-Z operators. The field-offset vector contains 512 points from -0.08 to +0.08; the plot labels this offset axis in Tesla. The powder grid is `rep_2ang_6400pts_sph`, the steady-state solution method is `backslash`, and the sequence requests `aniso_eq`.

The calculation is `powder(spin_system,@dnp_field_scan,parameters,'esr')`. The resulting spectrum is made real for plotting against `parameters.fields`, with the 15N S_z expectation value as the ordinate.

## Result and scope

The no-argument function creates a field-profile plot and does not declare a MATLAB return value. Its output is tied to the specified two-spin model, relaxation settings, single electron drive setting, 512-point field grid, and spherical powder grid; the source gives no numerical spectrum in the file. It depends on Spinach and its `dnp_field_scan`, steady-state, and powder-averaging routines.
