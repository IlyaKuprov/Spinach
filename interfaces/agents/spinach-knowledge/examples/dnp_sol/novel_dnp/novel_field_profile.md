# examples/dnp_sol/novel_dnp/novel_field_profile.m

- Signature: `novel_field_profile()`
- Source: [`examples/dnp_sol/novel_dnp/novel_field_profile.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/novel_dnp/novel_field_profile.m)

## Purpose

Computes the proton longitudinal signal at the end of a NOVEL contact sequence over a one-dimensional microwave-offset scan. The source describes the observable as the 1H I_z expectation value after 0.25 μs and cites [Redrouthu et al., DOI 10.1063/1.5000528](https://doi.org/10.1063/1.5000528). The source estimates seconds for the calculation.

## Spin model

The function sets `sys.magnet=0.34` (described in the source as an X-band magnet) and models one electron and two 1H spins. The electron principal g values are `[2.00319 2.00319 2.00258]`; the two proton Zeeman entries are the source's “ppm guess” tensors `[0 0 5]` and `[0 5 0]`. Their Euler-angle entries are supplied as `(pi/180)*{[0 10 0],[0 0 10],[100 0 0]}`. The three Cartesian coordinate rows are `[0 0 0]`, `[0 3.5 0]`, and `[2.475 2.475 0]`; the source does not label their units. It sets `inter.temperature=80` without annotating a unit.

The basis is `zeeman-hilb` with `approximation='none'`. Detection uses the proton state `state(spin_system,'Lz','1H')`.

## Sequence and offset scan

The NOVEL parameters select spins `{'E','1H'}`, `flippulse=1`, a timestep of `1e-9` seconds, and 250 steps (the source's 0.25 μs contact time). The electron nutation-frequency parameter is `14.48e6`; the 90-degree pulse duration is set by `1/(4*parameters.irr_powers)`. The powder grid is `rep_2ang_100pts_sph`, and the sequence requests `aniso_eq`.

The input offset array is 71 points from -35e6 to +35e6 Hz. Each point is shifted by `reference_point=-3.3e6` before being passed as `localpar.offset=[offset 0]`. For each offset the function calls `powder(spin_system,@noveldnp,localpar,'esr')` and stores the real part of the final returned signal element. These offset simulations run in `parfor`.

## Result and scope

The function plots the stored final proton signal against the unshifted offset axis converted to MHz; the ordinate is labelled as the 1H I_z expectation value. It returns no MATLAB output argument: the result is the figure. This is a fixed three-spin, powder-averaged scan on the stated grid and offset range, rather than a parameterised driver for arbitrary systems or pulse settings. It depends on Spinach, its `noveldnp` sequence and powder machinery, and MATLAB support for `parfor`.
