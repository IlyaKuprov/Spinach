# examples/esr_sol_pulsed/hard_3_pulse_deer_gd_1.m

- Signature: `hard_3_pulse_deer_gd_1()`

## Purpose

Gadolinium(III) DEER experiment at W-band using ideal pulses. Set to reproduce Figure 2b from Otting and co-authors: http://dx.doi.org/10.1021/ja204415w. The calculation uses brute-force time propagation and grid powder averaging, with central transitions on both gadolinium ions. Calculation time: minutes.

## Physical / mathematical content

- The system contains two `E8` spins at [0.00, 0.00, 0.00] and [60.50, 0.00, 0.00], with each zero-field-splitting matrix diag(1e8, 1e8, -2e8), at a 3.5 T magnetic field; the isotropic electron g-factor is 2.002319.
- The source notes that simulated gadolinium spin echoes are difficult to catch without the zero-field-splitting distributions found experimentally. Flip-flop terms in the inter-electron dipolar interaction are switched off using `deer-zz` to mimic slightly different experimental pulse frequencies.

## Numerical / algorithmic content

- The sequence uses a 100 ns step, 80 steps, and the `rep_2ang_1600pts_sph` powder grid. Its spectrum settings are a 1e10 sweep parameter and 1024 steps.
- Three pulse FIDs receive exponential apodisation with parameter 6 before FFTs using four times the spectrum step count.

## Implementation structure

- Create the spin system in the `zeeman-hilb` basis without approximation; run `powder` with `@deer_3p_hard_deer` in the `deer-zz` context.
- Plot the frequency-swept spectrum, probe and pump excitation profiles, and the negative imaginary DEER trace.
