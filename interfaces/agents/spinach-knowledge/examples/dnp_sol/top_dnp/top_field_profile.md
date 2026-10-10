# examples/dnp_sol/top_dnp/top_field_profile.m

- MATLAB implementation: [examples/dnp_sol/top_dnp/top_field_profile.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/top_dnp/top_field_profile.m)

## Purpose

Calculates a TOP-DNP field profile: the proton `I_z` expectation value at the end of a fixed contact time as microwave resonance offset is varied. The source estimates minutes of calculation time for a large powder grid. Background: [Redrouthu et al., *Science Advances* (2019)](https://doi.org/10.1126/sciadv.aav6909).

## Run and model

Call `top_field_profile()` in MATLAB. The source uses a Q-band setting `sys.magnet=1.2142`, one electron and two protons, spin temperature 80, and `zeeman-hilb` basis with no approximation. The trityl electron Zeeman values are `[2.00319 2.00319 2.00258]`; the proton values are `[0 0 5]` and `[0 5 0]`; the source describes the trityl values as a g-tensor and the 1H values as ppm guesses. Euler angles are `[0 10 0]`, `[0 0 10]`, and `[100 0 0]` degrees. The coordinates are `[0 0 0]`, `[0 3.5 0]`, and `[2.475 2.475 0]`; their units are not specified in the source.

The proton `L_z` state is detected. Each point uses `powder(spin_system,@topdnp,localpar,'esr')` with `parameters.needs={'aniso_eq'}` (source comment: the sequence needs `rho_eq`) and the spherical grid `rep_2ang_3200pts_sph`.

## Dependencies

Uses Spinach `create`, `basis`, `state`, and `powder` with `topdnp`; the offset loop is written as MATLAB `parfor`. Plotting uses `kfigure`, `plot`, and Spinach axis helpers.

## Offset profile

The electron nutation frequency is held at `17.8e6` Hz; the pulse and delay are 10 ns and 14 ns, respectively, with 300 TOP-DNP blocks. The 120-point offset grid is `linspace(-150e6,150e6,120)`, with `reference_point=-13e6` added to the first offset component for each run. Thus the varying scan dimension is resonance offset; the source does not sweep pulse amplitude in this function.

## Output and scope

A `parfor` loop evaluates the contact curve at each offset and stores `real(contact_curve(end))` as the profile value. The line plot uses the 120 scan offsets in MHz on the horizontal axis and the proton `I_z` expectation value on the vertical axis. This profile is limited to the configured offset grid, fixed pulse settings, TOP sequence, and powder grid.
