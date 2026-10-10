# examples/dnp_sol/beam_dnp/beam_parameter_scan.m

- MATLAB implementation: [examples/dnp_sol/beam_dnp/beam_parameter_scan.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/beam_dnp/beam_parameter_scan.m)

- **Call:** `beam_parameter_scan()` (zero input arguments; the function updates a contour figure and returns no explicit value).

## Purpose

Maps the final proton $I_z$ expectation from a BEAM DNP contact over microwave resonance offset and electron pulse nutation frequency. The cited experimental context is [Redrouthu et al., *Science Advances*](https://doi.org/10.1126/sciadv.abq0536); the source estimates minutes for its large powder-grid calculation.

## Model and source parameters

The X-band model sets `sys.magnet=0.3483` and contains one electron plus two protons (`{'E','1H','1H'}`). It uses the trityl g-tensor `[2.00319 2.00319 2.00258]`, proton entries `[0 0 5]` and `[0 5 0]` (source-described ppm guesses), Euler entries `[0 10 0]`, `[0 0 10]`, and `[100 0 0]` degrees, and Cartesian coordinates `[0,0,0]`, `[0,3.500,0]`, `[2.475,2.475,0]`. The coordinate and `inter.temperature=80` units are not stated in the source. The full `zeeman-hilb` basis has no approximation. Detection is `coil_state(spin_system,'Lz','1H','exact')`.

The sequence uses `parameters.spins={'E','1H'}`, pulse durations `20.0e-9` and `28.7e-9` seconds, `165` BEAM blocks, the `rep_2ang_800pts_sph` powder grid, and `parameters.needs={'aniso_eq'}` (source comment: “Sequence needs rho_eq”).

## Scan and output

The offset axis has `120` points from `-60e6` to `60e6` Hz; `reference_point=-3.3e6` is added to each offset. The electron nutation-frequency axis has `30` points from `20e6` to `40e6` Hz. For each offset, a `parfor` loop sets `irr_powers` to a nutation-frequency value, calls `powder(spin_system,@beamdnp,localpar,'esr')`, and stores the real part of the final contact-curve point. The plotted axes are offset and nutation frequency in MHz; the contour surface has 100 levels and is refreshed as each offset column is calculated.

**Dependencies:** Spinach system, basis, state, powder and plotting routines; the `beamdnp` sequence and named powder grid; MATLAB parallel-loop support for `parfor`.

**Use when:** the dependence on both microwave offset and electron nutation frequency is required. The script does not return the calculated surface as a function output; it builds the local `dnp_surf` array and presents it in the contour plot.
