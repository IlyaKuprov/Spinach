# examples/dnp_sol/novel_dnp/novel_contact_curve.m

- Signature: `novel_contact_curve()`

## Purpose

Calculates the contact-time evolution for the NOVEL solid-effect DNP example, described in the source as the transformation of (-E_z) into (I_z). It plots the real proton (I_z) expectation value against contact time. Further information is cited at [doi:10.1063/1.5000528](https://doi.org/10.1063/1.5000528). The source estimates the calculation time as seconds.

## Physical / mathematical content

- The model has one electron and two protons at 0.34 T and 80 K. The electron g-tensor principal values are 2.00319, 2.00319, and 2.00258; the proton Zeeman values are zero with the source's stated ppm guesses.
- The listed positions are ([0,0,0]), ([0,3.5,0]), and ([2.475,2.475,0]). The source uses a full Zeeman-Hilbert basis.
- The NOVEL calculation uses an electron offset of -3.3 MHz, a 14.48 MHz electron pulse nutation frequency, 1 ns time steps, 2,400 steps, and the `rep_2ang_400pts_sph` powder grid.

## Numerical / algorithmic content

Creates and bases the Spinach system with `formalism='zeeman-hilb'` and `approximation='none'`; calls `powder` with the `noveldnp` callback in ESR mode; and plots the real part of the returned curve over the contact-time axis. The sequence is configured as a flip pulse and requests the anisotropic-equilibrium term.

## Implementation structure

The function specifies the field, isotopes, anisotropic Zeeman parameters and orientations, coordinates, and temperature; builds the system and basis; sets proton detection and the NOVEL pulse, offset, time-step, and powder-grid parameters; runs the powder calculation; and plots the proton (I_z) signal against time.
