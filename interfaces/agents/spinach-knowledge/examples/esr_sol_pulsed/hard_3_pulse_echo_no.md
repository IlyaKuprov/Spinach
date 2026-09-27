# examples/esr_sol_pulsed/hard_3_pulse_echo_no.m

- Signature: `hard_3_pulse_echo_no()`

## Purpose

Simulates a three-pulse DEER spin echo for two nitroxide radicals separated by 25 Å at X-band, using brute-force time propagation and numerical powder averaging.

## Model and calculation

The electron spins are simulated at 0.33 T with anisotropic nitroxide g tensors and the orientations given in the source. The calculation uses a Zeeman Hilbert-space basis and a 3,200-point spherical powder grid—not a Liouville-space basis. It evaluates the hard-echo sequence with delays of 1.0 µs, 0.5 µs, and 0.5 µs, using 256 steps, then plots the echo against time in microseconds. The source estimates a run time of seconds.

Nitroxide g-tensor reference: http://dx.doi.org/10.1063/1.1697233
