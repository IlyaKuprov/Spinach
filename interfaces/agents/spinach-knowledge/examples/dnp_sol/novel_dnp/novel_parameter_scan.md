# examples/dnp_sol/novel_dnp/novel_parameter_scan.m

- Signature: `novel_parameter_scan()`

## Purpose

Computes a two-dimensional NOVEL DNP parameter scan: the proton (I_z) expectation value at the end of a 0.25 (mumathrm{s}) contact sequence as a function of electron-pulse nutation frequency and microwave resonance offset. The source cites [doi:10.1063/1.5000528](https://doi.org/10.1063/1.5000528) and estimates the calculation time as minutes.

## Physical / mathematical content

- The system is one electron and two protons at 0.34 T and 80 K, with the trityl-like electron g-tensor, proton Zeeman guesses, and coordinates specified in the companion NOVEL examples.
- The calculation uses 250 steps of 1 ns, a `rep_2ang_100pts_sph` powder grid, and anisotropic equilibrium. For each nutation frequency it sets the pulse duration to a quarter cycle.
- The scan covers 30 electron nutation frequencies from 1 to 30 MHz and 71 offsets from -35 to +35 MHz, shifted by the -3.3 MHz reference point. The output surface stores the real final-time proton (I_z) expectation value.

## Numerical / algorithmic content

Uses a full Zeeman-Hilbert basis and runs the `noveldnp` sequence through `powder` in ESR mode. A serial outer loop sets each nutation frequency and pulse duration; a `parfor` loop evaluates its microwave offsets. The result is updated as a 100-level filled contour plot with offset and nutation frequency axes.

## Implementation structure

The function builds the electron/two-proton system and basis, configures proton detection and the NOVEL sequence, generates the offset and nutation-frequency axes, then nests the per-power setup around the parallel offset simulations. Each simulation contributes its final proton signal to the surface before the contour plot is refreshed.
