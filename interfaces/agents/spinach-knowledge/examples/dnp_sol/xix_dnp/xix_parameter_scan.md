# examples/dnp_sol/xix_dnp/xix_parameter_scan.m

- MATLAB implementation: [examples/dnp_sol/xix_dnp/xix_parameter_scan.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/xix_dnp/xix_parameter_scan.m)

## Purpose

`xix_parameter_scan()` maps the final proton longitudinal-polarisation signal of a XiX DNP contact over microwave resonance offset and electron nutation frequency. The source links the example to [the XiX DNP study](https://doi.org/10.1021/jacs.1c09900).

## Model and sequence

The model is a trityl electron and two protons at a Q-band setting (`sys.magnet=1.2142`), with trityl g principal values `[2.00319 2.00319 2.00258]` and proton Zeeman entries `[0 0 5]` and `[0 5 0]` (described in the source as ppm guesses). Its orientation triples are `[0 10 0]`, `[0 0 10]`, and `[100 0 0]` degrees, converted to radians in the script. The Cartesian coordinates are `[0 0 0]`, `[0 3.5 0]`, and `[2.475 2.475 0]`; their units are not stated. Spin temperature is set to `80` without an explicit unit. The basis is `zeeman-hilb` with no approximation, and detection uses proton `Lz`.

The sequence is `@xixdnp` evaluated through `powder(...,'esr')`. It uses `{'E','1H'}`, 150 XiX blocks, `48e-9`-second pulses, phase `pi` (the second pulse is described as opposite phase), spherical grid `rep_2ang_400pts_sph`, and `needs={'aniso_eq'}` (the source comment says the sequence needs `rho_eq`). Nutation frequency is assigned through `irr_powers`.

## Two-parameter scan

The scan contains 120 nominal offsets from -100e6 to +100e6 Hz and 30 electron nutation frequencies from 10e6 to 50e6 Hz. A -13e6 Hz reference is added to each simulated first offset; the second offset component is zero. The offset axis is displayed in MHz. For every grid point the plotted quantity is the real final contact-curve value, rather than the full time-dependent curve. The script loops over offsets and uses `parfor` over nutation frequencies within each offset.

## Use and output

With Spinach available on the MATLAB path, call `xix_parameter_scan()`. It evaluates the powder-averaged XiX contact over both parameter grids and displays a 100-level contour plot with offset and electron nutation frequency in MHz. The function declares no return value and does not save a scan matrix; its result is the figure. This scan changes offset and nutation frequency only, with the model, sequence length, pulse duration, phase, and powder grid held fixed.

The example depends on Spinach system/basis/state construction, the `xixdnp` sequence, `powder` ESR averaging, and Spinach plotting helpers; its inner scan uses MATLAB `parfor`.
