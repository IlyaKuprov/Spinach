# examples/dnp_sol/beam_dnp/beam_field_profile.m

- MATLAB implementation: [examples/dnp_sol/beam_dnp/beam_field_profile.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/beam_dnp/beam_field_profile.m)

- **Call:** `beam_field_profile()` (zero input arguments; the function produces a figure rather than a returned value).

## Purpose

Calculates the proton $I_z$ expectation at the end of a BEAM DNP contact as a function of microwave resonance offset, with electron pulse amplitude held at the source-set value. The example cites [Redrouthu et al., *Science Advances*](https://doi.org/10.1126/sciadv.abq0536). Its source estimates minutes of calculation because it uses a large powder grid.

## Model and source parameters

The X-band model sets `sys.magnet=0.3483`, with one electron and two protons (`{'E','1H','1H'}`). The electron trityl g-tensor is `[2.00319 2.00319 2.00258]`; the two proton Zeeman entries are `[0 0 5]` and `[0 5 0]`, described in the source as ppm guesses. The Euler-angle entries are `[0 10 0]`, `[0 0 10]`, and `[100 0 0]` degrees. Cartesian coordinates are `[0,0,0]`, `[0,3.500,0]`, and `[2.475,2.475,0]`; the source gives no coordinate unit. The source sets `inter.temperature=80` without stating a unit.

The full Zeeman-Hilbert basis uses `bas.formalism='zeeman-hilb'` and `bas.approximation='none'`. The detected state is `state(spin_system,'Lz','1H')`. BEAM sequence parameters are `parameters.spins={'E','1H'}`, electron nutation frequency `32.0 MHz` (documented as `32.0e6 Hz`), pulse durations `20.0 ns` and `28.7 ns`, and `165` BEAM blocks. Powder averaging uses `rep_2ang_800pts_sph`; `parameters.needs={'aniso_eq'}` is accompanied by the source comment “Sequence needs rho_eq”.

## Scan, calculation, and output

The script forms `120` offsets from `-60e6` to `60e6` Hz and adds `reference_point=-3.3e6` to each electron offset. A `parfor` loop calls `powder(spin_system,@beamdnp,localpar,'esr')` at each offset; it stores the real part of the final contact-curve value. The plot converts the scan axis to MHz and labels the vertical observable as the $I_z$ expectation on $^1$H.

**Dependencies:** Spinach system/basis/state and powder routines, the `beamdnp` sequence and named powder grid, MATLAB parallel-loop support for `parfor`, and the Spinach plotting helpers `kfigure`, `kylabel`, `kxlabel`, and `kgrid`.

**Use when:** an offset profile at the specified fixed BEAM pulse amplitude is wanted. The source does not expose the profile array as a function return value or set a separately named contact-time parameter; the endpoint is the last value returned by the sequence calculation.
