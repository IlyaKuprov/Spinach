# examples/dnp_sol/tppm_dnp/tppm_contact_curve.m

- MATLAB implementation: [examples/dnp_sol/tppm_dnp/tppm_contact_curve.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/tppm_dnp/tppm_contact_curve.m)

**Entry point:** tppm_contact_curve() (no input arguments).

## Purpose

Calculates and plots the proton z-polarisation expectation during a TPPM DNP contact, following transfer from electron −E_z to proton I_z. The source lists a calculation time of seconds.

## Spin system and sequence

The example sets a Q-band magnet with **sys.magnet=1.2142** and a three-spin system comprising one electron and two protons. The electron Zeeman principal values are [2.00319, 2.00319, 2.00258]; the proton Zeeman values are [0, 0, 5] and [0, 5, 0] (the source labels these proton values ppm guesses). The Euler-angle arrays are [0, 10, 0], [0, 0, 10], and [100, 0, 0] degrees, converted to radians in the script. Cartesian coordinates are [0,0,0], [0,3.5,0], and [2.475,2.475,0]. The spin temperature is set to 80. The basis is **zeeman-hilb** with **approximation='none'**.

The experiment uses **spins={'E','1H'}**, an electron nutation frequency of 33 MHz (**irr_powers=33e6** Hz), 16 ns pulses, 250 TPPM blocks, and a second-pulse phase of 120°. The sequence offset is [−9e6, 0] Hz, formed from the −13 MHz reference point plus a 4 MHz offset. Powder averaging uses the 400-point **rep_2ang_400pts_sph** grid.

## Calculation and output

After creating the spin system and basis, the script detects proton **Lz** and calls **powder(spin_system,@xixdnp,parameters,'esr')**; **parameters.needs** is **{'aniso_eq'}**. The returned contact curve is plotted as its real part against a time vector from zero to **2*pulse_dur*nloops** seconds (8 μs, with **nloops+1** points). The figure labels the ordinate as the proton I_z expectation value. The script displays a plot rather than writing a data file.

The calculation depends on Spinach system/basis/state construction, the **powder** averaging routine, the **xixdnp** sequence, and the **kfigure**, **kxlabel**, and **kylabel** plotting helpers.

## Reference

[Redrouthu et al., DOI: 10.1063/5.0153053](https://doi.org/10.1063/5.0153053).
