# examples/esr_sol_pulsed/hard_3_pulse_deer_exchange.m

- Signature: `hard_3_pulse_deer_exchange()`

## Purpose

Three-pulse DEER on a Cu(II)–Cu(II) system in a linked porphyrin complex with strong exchange coupling between the electrons. A distribution of exchange couplings is summed over. The calculation uses brute-force time propagation and numerical powder averaging in Liouville space. Calculation time: minutes.

## Physical / mathematical content

- Two electrons with g eigenvalues [2.050, 2.050, 2.195] are at [0, 0, 0] and [24.50, 0, 0] in a 1.2132 T magnetic field.
- Exchange couplings span 6e6 to 20e6 in 20 steps; Gaussian weights centred at 13.1e6 with parameter 4.2622e6 are normalised before summation.

## Numerical / algorithmic content

- Each simulation uses a 2.5 ns step, 200 steps, and the `rep_2ang_400pts_sph` powder grid.

## Implementation structure

- For each exchange coupling, create the spin system in the `sphten-liouv` basis without approximation, then call `powder` with `@deer_3p_hard_deer` in the `deer-zz` context.
- Sum the weighted DEER traces and plot the imaginary part against time in microseconds.
