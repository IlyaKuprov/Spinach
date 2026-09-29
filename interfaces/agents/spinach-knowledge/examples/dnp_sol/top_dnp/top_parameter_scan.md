# examples/dnp_sol/top_dnp/top_parameter_scan.m

- MATLAB implementation: [examples/dnp_sol/top_dnp/top_parameter_scan.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/top_dnp/top_parameter_scan.m)

## Purpose

Maps the terminal proton `I_z` expectation value for a TOP-DNP experiment against microwave resonance offset and electron pulse nutation frequency, at fixed pulse and delay durations. The source estimates minutes of calculation time for a large powder grid. Background: [Redrouthu et al., *Science Advances* (2019)](https://doi.org/10.1126/sciadv.aav6909).

## Run and model

Call `top_parameter_scan()` in MATLAB. The model uses a Q-band setting `sys.magnet=1.2142`, one electron and two protons, spin temperature 80, and a `zeeman-hilb` basis without approximation. The trityl electron Zeeman values are `[2.00319 2.00319 2.00258]`; proton values are `[0 0 5]` and `[0 5 0]`; the source describes the trityl values as a g-tensor and the 1H values as ppm guesses. Euler angles `[0 10 0]`, `[0 0 10]`, and `[100 0 0]` degrees. The coordinates are `[0 0 0]`, `[0 3.5 0]`, and `[2.475 2.475 0]` (units are not specified in the source).

The calculation detects proton `L_z` and calls `powder(spin_system,@topdnp,localpar,'esr')` for each scan point. It uses `parameters.needs={'aniso_eq'}` (annotated as requiring `rho_eq`) and spherical grid `rep_2ang_400pts_sph`.

## Dependencies

Uses Spinach `create`, `basis`, `state`, and `powder` with `topdnp`; plotting uses `kfigure`, `contourf`, and Spinach axis/colorbar helpers.

## Scan settings

The two axes are 120 resonance offsets from `-100e6` to `100e6` and 30 electron nutation frequencies from `10e6` to `50e6`. The source labels these grids in Hz and plots both axes in MHz; `reference_point=-13e6` is added to the first offset component. Pulse duration is 10 ns, delay duration is 14 ns, and `nloops=300`. For each pair, the code sets `localpar.irr_powers` to the selected nutation frequency and records the real part of the final contact-curve point.

## Output and scope

The values fill a `numel(nutfrqs)` by `numel(offsets)` surface and are drawn as a 100-level filled contour: offset (MHz) horizontally, electron nutation frequency (MHz) vertically, and proton `I_z` expectation value by colour. The script plots the finite 120-by-30 grid; it does not save a separate numeric result file. Its TOP sequence settings and orientation grid remain fixed across that scan.
