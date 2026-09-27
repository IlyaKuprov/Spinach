# examples/esr_sol_pulsed/hard_3_pulse_echo_gd.m

- Signature: `hard_3_pulse_echo_gd()`

## Purpose

Calculates a three-pulse DEER echo for two Gd(III) centres. The outermost ZFS transition is excited by the probe pulse and the central transition by the pump pulse; the pulses are treated as hard. The calculation uses brute-force time propagation and powder averaging.

## Model and calculation

The pair is modelled at 3.5 T with axial ZFS tensors of magnitude `D = 0.56 GHz`, different tensor orientations, and a 29.5 Å separation. The pulse delays are 2 µs, 1 µs, and 50 ns; the trace has 500 steps and is averaged over a 1,600-point spherical grid. The script plots the imaginary echo as a function of time in microseconds.

The source cautions that the Gd echo is very sharp and that the simulation does not include the ZFS distributions found in experimental systems. Estimated calculation time: seconds.
