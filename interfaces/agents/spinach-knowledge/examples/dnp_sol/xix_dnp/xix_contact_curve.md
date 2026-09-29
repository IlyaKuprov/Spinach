# examples/dnp_sol/xix_dnp/xix_contact_curve.m

- MATLAB implementation: [examples/dnp_sol/xix_dnp/xix_contact_curve.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/xix_dnp/xix_contact_curve.m)

**Entry point:** xix_contact_curve() (no input arguments).

## Purpose

Calculates and plots the proton I_z expectation during an X-inverse-X (XiX) DNP contact, following transfer from electron −E_z to proton I_z. The source lists a calculation time of seconds.

## Spin system and sequence

The Q-band model contains one electron and two protons, with **sys.magnet=1.2142**, electron Zeeman principal values [2.00319, 2.00319, 2.00258], proton Zeeman values [0, 0, 5] and [0, 5, 0] (identified as ppm guesses in the source), and Cartesian coordinates [0,0,0], [0,3.5,0], and [2.475,2.475,0]. The spin temperature is set to 80, with a full **zeeman-hilb** basis and no approximation.

The experiment uses **spins={'E','1H'}**, electron nutation frequency 17.8 MHz (**irr_powers=17.8e6** Hz), 48 ns pulses, 80 XiX blocks, an inverted second-pulse phase of π, and the 1600-point **rep_2ang_1600pts_sph** powder grid. Its offset argument is [46.5e6,0] Hz: the source forms this from the −13 MHz reference point and a 59.5 MHz offset. The script sets **parameters.needs={'aniso_eq'}**.

## Calculation and output

The script constructs the spin system and basis, detects proton **Lz**, then calls **powder(spin_system,@xixdnp,parameters,'esr')**. It plots the real part of the returned contact curve against a time vector from zero to **2*pulse_dur*nloops** seconds (7.68 μs, with **nloops+1** points) and labels the ordinate as proton I_z expectation.

Dependencies are Spinach system, basis, and state construction; **powder** and the **xixdnp** sequence; and the **kfigure**, **kxlabel**, and **kylabel** plotting helpers.

## Reference

[XiX DNP paper, DOI: 10.1021/jacs.1c09900](https://doi.org/10.1021/jacs.1c09900).
