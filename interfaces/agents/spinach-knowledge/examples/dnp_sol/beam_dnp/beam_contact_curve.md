# examples/dnp_sol/beam_dnp/beam_contact_curve.m

- MATLAB implementation: [examples/dnp_sol/beam_dnp/beam_contact_curve.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/beam_dnp/beam_contact_curve.m)

## Purpose

Run `beam_contact_curve()` to plot the proton-`I_z` signal over the BEAM DNP contact sequence, illustrating the source-described transformation of `-E_z` into `I_z`. Further information is in the [Science Advances paper](https://doi.org/10.1126/sciadv.abq0536); the source header estimates seconds.

## Spin system and sequence

The X-band field is assigned `0.3483`. The spins are one electron and two `^1H` nuclei. The electron g eigenvalues are `[2.00319 2.00319 2.00258]`; proton Zeeman eigenvalues are `[0 0 5]` and `[0 5 0]`. The Euler-angle entries are `[0 10 0]`, `[0 0 10]`, and `[100 0 0]`, multiplied by `pi/180`. Coordinates are `[0 0 0]`, `[0 3.500 0]`, and `[2.475 2.475 0]`; temperature is assigned `80`. The basis is full Zeeman-Hilbert (`formalism='zeeman-hilb'`, `approximation='none'`), and detection is `coil_state(spin_system,'Lz','1H','exact')`.

The sequence parameters are spins `{'E','1H'}`, offset assignment `[(-3.3+5.0)*1e6 0]`, electron nutation frequency `32.0e6 Hz`, pulse durations `[20.0e-9 28.7e-9]` seconds, `nloops=165`, powder grid `'rep_2ang_800pts_sph'`, and `needs={'aniso_eq'}` (the source comment says the sequence needs `rho_eq`). The offset assignment's adjacent comment reads “-13 MHz reference point, 5.0 MHz offset”; the first-component expression itself evaluates to `+1.7e6`, so the source leaves a discrepancy between expression and comment.

The powder-averaged simulation call is `contact_curve=powder(spin_system,@beamdnp,parameters,'esr')`. The plotted time samples span zero to `sum(parameters.pulse_dur)*parameters.nloops` using `nloops+1` points; the plot is `real(contact_curve)`, with contact time in seconds and the `^1H I_z` expectation on the vertical axis.
