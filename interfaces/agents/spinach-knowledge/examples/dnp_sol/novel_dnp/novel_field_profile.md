# examples/dnp_sol/novel_dnp/novel_field_profile.m

- Signature: `novel_field_profile()`

## Purpose

Calculates the field profile of a NOVEL DNP experiment: the proton (I_z) expectation value after a 0.25 μs contact period, as a function of electron-pulse resonance offset. The source cites [doi:10.1063/1.5000528](https://doi.org/10.1063/1.5000528) and estimates the calculation time as seconds.

## Physical / mathematical content

- The system contains one electron and two protons at 0.34 T and 80 K, with the same trityl-like electron g-tensor and proton Zeeman guesses as the companion NOVEL examples. The source specifies coordinates ([0,0,0]), ([0,3.5,0]), and ([2.475,2.475,0]).
- Each simulated contact sequence has 250 steps of 1 ns, a fixed electron nutation frequency of 14.48 MHz, and a `rep_2ang_100pts_sph` powder grid.
- The scan evaluates 71 offsets from -35 to +35 MHz, added to a -3.3 MHz reference point. At each offset the plotted datum is the final proton (I_z) expectation value.

## Numerical / algorithmic content

Uses a full Zeeman-Hilbert basis and the `noveldnp` sequence through `powder` in ESR mode. A `parfor` loop independently simulates each microwave offset; the script then plots the real final-time signal against offset in MHz. Output is hushed during the calculation.

## Implementation structure

The function builds the electron/two-proton system and basis, configures proton detection and the NOVEL pulse, time-step, and powder-grid parameters, loops over the frequency offsets while setting each local parameter structure's offset, extracts the final signal point, and plots the field profile.
