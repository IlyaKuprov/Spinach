# examples/dnp_sol/tppm_dnp/tppm_field_profile.m

- MATLAB implementation: [examples/dnp_sol/tppm_dnp/tppm_field_profile.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/tppm_dnp/tppm_field_profile.m)

**Entry point:** tppm_field_profile() (no input arguments).

## Purpose

Calculates the final proton I_z expectation after a fixed-duration TPPM DNP contact as a function of microwave resonance offset. Here “field profile” is the frequency-offset profile shown on the x-axis, not a sweep of the static magnet setting. The source estimates minutes for the powder-averaged calculation.

## Spin system and sequence

The model is the Q-band three-spin system (one electron and two protons), with **sys.magnet=1.2142**, electron Zeeman principal values [2.00319, 2.00319, 2.00258], proton Zeeman values [0, 0, 5] and [0, 5, 0] (described in the source as ppm guesses), and Cartesian coordinates [0,0,0], [0,3.5,0], and [2.475,2.475,0]. Its spin temperature is set to 80; it uses the full **zeeman-hilb** basis with no approximation. Output from system setup is hushed.

Each profile point uses an electron nutation frequency of 17.8 MHz (**irr_powers=17.8e6** Hz), 48 ns pulses, 150 TPPM blocks, second-pulse phase 120°, and the 1600-point **rep_2ang_1600pts_sph** spherical grid. The sequence uses **spins={'E','1H'}** and **parameters.needs={'aniso_eq'}**. A 120-point offset grid spans −150 to 150 MHz; each sequence offset adds the −13 MHz reference point to the swept value.

## Calculation and output

The script parallelises over offset values with **parfor**. At each point it calls **powder(spin_system,@xixdnp,localpar,'esr')** and stores **real(contact_curve(end))**, the final value of the returned contact curve. It plots those final values against the swept offsets converted to MHz and labels the ordinate as proton I_z expectation. With 150 blocks of 48 ns pulses, the plotted value is taken at the end of a 14.4 μs contact.

Dependencies are Spinach system, basis, and state construction; **powder** and the **xixdnp** sequence; MATLAB **parfor**; and the **kfigure**, **kxlabel**, and **kylabel** plotting helpers.

## Reference

[Redrouthu et al., DOI: 10.1063/5.0153053](https://doi.org/10.1063/5.0153053).
