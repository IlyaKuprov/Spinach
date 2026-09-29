# examples/optimal_control/pattern_pulse_1.m

[Source file](https://github.com/IlyaKuprov/Spinach/blob/main/examples/optimal_control/pattern_pulse_1.m)

## Purpose

This example designs a phase-modulated, nutation-frequency-selective excitation pulse. For an on-resonance 13C spin at 28.18 T with no basis approximation, it asks for an initial Sz state to reach Sz in three selected nutation-frequency bands and Sx at the other sampled frequencies. The source cites the Glaser-group paper at [doi:10.1016/j.jmr.2004.12.005](https://doi.org/10.1016/j.jmr.2004.12.005).

## Target pattern and pulse design

The drift Hamiltonian is zero for this single on-resonance spin. The target pattern is defined over 128 nutation frequencies from 6 to 14 kHz: Sz is targeted at indices 1-20, 54-73, and 109-128; Sx is targeted at the remaining samples. The source plots this requested pattern before optimisation. It keeps the amplitude profile fixed and optimises phase with GRAPE via `fmaxnewton` and `@grape_phase`, using the `lbfgs` method. The initial phase is constant at pi/4, the configured iteration limit is 200, and the pulse grid has 250 intervals of 20 microseconds (5 milliseconds total). The B1 levels are set from the 6-14 kHz range and represented in the controls as angular frequencies in rad/s.

## Evaluation shown by the example

For each nutation-frequency sample, the script simulates the shaped pulse and projects the final state onto Sx and Sz, then plots those projections against the target pattern. This defines how the example evaluates the design; the source contains no measured or reported post-optimisation values. Its calculation-time comment says minutes. Source contact: ilya.kuprov@weizmann.ac.il.
