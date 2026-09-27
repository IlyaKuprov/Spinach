# examples/esr_sol_pulsed/hard_3_pulse_deer_cu.m

- Signature: `hard_3_pulse_deer_cu()`

## Purpose

Three-pulse DEER on a Cu(II)–NO two-electron system at X-band. The numerical calculation uses brute-force time propagation and powder averaging in Liouville space, including orientation effects of the electron g tensors on dipolar coupling. An analytical calculation uses the isotropic parts of the electron g tensors. Calculation time: seconds.

## Physical / mathematical content

- The electrons are at [0, 0, 0] and [20, 0, 0] in a 0.33 T magnetic field, with g eigenvalues [2.056, 2.056, 2.205] and [2.009, 2.006, 2.003].
- The analytical dipolar coupling is calculated from the coordinates and scaled using the mean g-factor of each electron.

## Numerical / algorithmic content

- The numerical sequence uses a 10 ns step, 50 steps, and the `rep_2ang_1600pts_sph` powder grid.
- The analytical trace is `0.35*deer_analyt(D,0,time_axis)`.

## Implementation structure

- Create the spin system in the `sphten-liouv` basis without approximation, with trajectory-level SSR disabled.
- Run `powder` with `@deer_3p_hard_deer` in the `deer` context, then plot the imaginary numerical trace beside the analytical result.
