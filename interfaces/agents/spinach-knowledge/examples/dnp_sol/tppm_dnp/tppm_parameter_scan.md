# examples/dnp_sol/tppm_dnp/tppm_parameter_scan.m

- MATLAB implementation: [examples/dnp_sol/tppm_dnp/tppm_parameter_scan.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/tppm_dnp/tppm_parameter_scan.m)

**Entry point:** tppm_parameter_scan() (no input arguments).

## Purpose

Maps the final proton I_z expectation after a fixed-duration TPPM DNP contact over electron nutation frequency and microwave resonance offset. The source estimates minutes for the powder-averaged scan.

## Spin system and sequence

The model is the Q-band system of one electron and two protons, with **sys.magnet=1.2142**, electron Zeeman principal values [2.00319, 2.00319, 2.00258], proton Zeeman values [0, 0, 5] and [0, 5, 0] (source-described ppm guesses), and coordinates [0,0,0], [0,3.5,0], and [2.475,2.475,0]. The spin temperature is set to 80. The script uses the **zeeman-hilb** basis with **approximation='none'**, and hushes setup output.

Every point uses 48 ns pulses, 150 TPPM blocks, a 120° second-pulse phase, the 400-point **rep_2ang_400pts_sph** powder grid, **spins={'E','1H'}**, and **parameters.needs={'aniso_eq'}**. The microwave offset grid has 120 points from −100 to 100 MHz; each sequence offset adds the −13 MHz reference point. The electron nutation-frequency grid has 30 points from 10 to 50 MHz.

## Calculation and output

For each offset, the script runs a **parfor** over nutation frequencies. Each task copies the experiment parameters, sets **irr_powers** in Hz, calls **powder(spin_system,@xixdnp,localpar,'esr')**, and records **real(contact_curve(end))**. Thus the surface contains the final proton signal for the selected contact duration (150 blocks × 2 × 48 ns = 14.4 μs), rather than the full contact curve. The updating contour plot displays microwave offset in MHz horizontally and electron nutation frequency in MHz vertically, with colour indicating proton I_z expectation.

Dependencies are Spinach system, basis, and state construction; **powder** and the **xixdnp** sequence; MATLAB **parfor**; and the **kfigure**, **contourf**, **kxlabel**, **kylabel**, and **kcolourbar** plotting helpers.

## Reference

[Redrouthu et al., DOI: 10.1063/5.0153053](https://doi.org/10.1063/5.0153053).
