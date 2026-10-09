# examples/dnp_sol/top_dnp/top_contact_curve.m

- MATLAB implementation: [examples/dnp_sol/top_dnp/top_contact_curve.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/top_dnp/top_contact_curve.m)

## Purpose

Plots the TOP-DNP contact curve: the source describes the transformation of `-E_z` into `I_z` during the contact time of the time-optimised pulsed DNP experiment. The source estimates seconds of calculation time. Background: [Redrouthu et al., *Science Advances* (2019)](https://doi.org/10.1126/sciadv.aav6909).

## Run and model

Call `top_contact_curve()` in MATLAB. The example builds a Q-band model with one electron and two protons, uses the trityl electron Zeeman values `[2.00319 2.00319 2.00258]` and proton values `[0 0 5]` and `[0 5 0]` (the source describes the trityl values as a g-tensor and the 1H values as ppm guesses), with Euler angles `[0 10 0]`, `[0 0 10]`, and `[100 0 0]` degrees. The coordinates are `[0 0 0]`, `[0 3.5 0]`, and `[2.475 2.475 0]` (the source does not label coordinate units); spin temperature is 80 and `sys.magnet=1.2142` is described as Q-band. The basis is `zeeman-hilb` with no approximation.

The detected state is the proton `L_z` state. The call `powder(spin_system,@topdnp,parameters,'esr')` evaluates the TOP sequence over the configured spherical grid `rep_2ang_3200pts_sph`; `parameters.needs={'aniso_eq'}` is annotated in the source as requiring `rho_eq`.

## Dependencies

Uses Spinach `create`, `basis`, `state`, and `powder` with the `topdnp` sequence callback; the plot uses `kfigure`, `plot`, Spinach axis helpers, and MATLAB plotting functions.

## Fixed sequence settings

The electron nutation frequency (`irr_powers`) is `17.8e6` Hz, pulse duration is 10 ns, delay duration is 14 ns, and the sequence has 300 TOP-DNP blocks. The microwave offset is assigned as `[(-13.0+92.5)*1e6 0]`; the source comment identifies a -13 MHz reference point and a 92.5 MHz offset. These remain fixed while this example evaluates the contact curve.

## Output and scope

The returned contact-curve values are plotted against a time axis from zero to `nloops*(pulse_dur+delay_dur)` seconds, with `nloops+1` points. The ordinate is the proton `I_z` expectation value. The script plots the curve and does not write a separate numeric data file. Its result is for this fixed TOP sequence configuration and powder grid.
