# examples/dnp_sol/novel_dnp/novel_contact_curve.m

- MATLAB implementation: [examples/dnp_sol/novel_dnp/novel_contact_curve.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/novel_dnp/novel_contact_curve.m)

## Purpose

This example follows the source's stated transformation of `-E_z` into `I_z` during contact in a NOVEL solid-effect DNP experiment. It calculates a contact-time curve and plots the real proton `I_z` expectation value. The cited reference is [doi:10.1063/1.5000528](https://doi.org/10.1063/1.5000528); the source estimates calculation time as seconds.

## Spin system and experiment

The no-argument function `novel_contact_curve()` sets an X-band magnet value of `0.34` and uses one electron and two protons. The electron g-tensor eigenvalues are `[2.00319,2.00319,2.00258]`; the proton shifts are source-labelled “ppm guess” values `[0,0,5]` and `[0,5,0]`. The Euler-angle lists are `[0,10,0]`, `[0,0,10]`, and `[100,0,0]` multiplied by `pi/180`. Coordinates are `[0,0,0]`, `[0,3.500,0]`, and `[2.475,2.475,0]`; coordinate units are not specified. The source sets `temperature=80` and uses the full `zeeman-hilb` basis (`approximation='none'`).

Detection is proton `Lz`. The experiment selects `spins={'E','1H'}`, sets `offset=[(-3.3+0.0)*1e6,0]` (commented as a -3.3 MHz reference point and 0.0 MHz offset), and uses an electron irradiation nutation frequency `irr_powers=14.48e6` Hz. The pulse duration is set to `1/(4*irr_powers)`; the sequence timestep is `1e-9` seconds, with 2,400 steps. It selects `rep_2ang_400pts_sph`, `flippulse=1`, and `needs={'aniso_eq'}`.

## Calculation and output

The calculation is `powder(spin_system,@noveldnp,parameters,'esr')`. The local `contact_curve` is plotted as `real(contact_curve)` versus a 2,401-point time axis from zero to `timestep*nsteps` (2.4 microseconds); the axis is labelled contact time in seconds. The function declares no output argument, so the displayed curve is the user-facing result. Running it requires Spinach and the `noveldnp` callback. The proton shifts are guesses in the source, and the result is tied to the specified finite spherical grid and pulse/offset setup.
