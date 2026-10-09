# examples/dnp_sol/xix_dnp/xix_field_profile.m

- MATLAB implementation: [examples/dnp_sol/xix_dnp/xix_field_profile.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/xix_dnp/xix_field_profile.m)

## Purpose

`xix_field_profile()` profiles the final proton longitudinal-polarisation signal from a XiX DNP contact against microwave resonance offset, with electron nutation frequency fixed. The source associates this example with [the XiX DNP study](https://doi.org/10.1021/jacs.1c09900).

## Model and sequence

The spin system is a trityl electron and two protons at a Q-band setting (`sys.magnet=1.2142`). The trityl g principal values are `[2.00319 2.00319 2.00258]`; the two proton Zeeman entries are `[0 0 5]` and `[0 5 0]`, which the source describes as ppm guesses. The orientation triples are `[0 10 0]`, `[0 0 10]`, and `[100 0 0]` degrees, converted to radians in the script. The Cartesian coordinates are E `[0 0 0]`, H `[0 3.5 0]`, and H `[2.475 2.475 0]` (coordinate units are not specified there). The spin-temperature setting is `80` (no unit is stated). The basis is `zeeman-hilb` with `approximation='none'`, and the detected operator is proton `Lz`.

The sequence is `@xixdnp` under `powder(...,'esr')`. Its settings select `{'E','1H'}`, electron nutation frequency `17.8e6` Hz, pulse duration `48e-9` seconds, 150 XiX blocks, phase `pi` (the source comment identifies the second pulse as opposite phase), spherical grid `rep_2ang_1600pts_sph`, and `needs={'aniso_eq'}` (the source comment says the sequence needs `rho_eq`).

## Offset profile

The nominal offset array has 120 points spanning -150e6 to +150e6 Hz. The simulated first offset is shifted by the -13e6 Hz reference point, with the second offset component set to zero; the plot labels the nominal, unshifted array in MHz. For each offset, the reported value is the real final point of the contact curve, not an average over its time points. Offsets are independent `parfor` work items.

## Use and output

With Spinach available on the MATLAB path, call `xix_field_profile()`. It builds the system and basis, evaluates the powder-averaged XiX contact at each offset, and opens a plot of the final proton `I_z` expectation value versus nominal microwave offset. The function declares no return value; the profile is used locally for plotting. This is a one-dimensional offset profile at fixed nutation frequency and fixed sequence/system settings, rather than a two-parameter optimisation or a returned dataset.

The calculation uses Spinach system/basis/state construction, the `xixdnp` sequence and `powder` ESR averaging, plus the Spinach plotting helpers. The offset loop uses MATLAB `parfor`.
