# examples/esr_sol_pulsed/hard_3_pulse_echo_cu.m

- Signature: `hard_3_pulse_echo_cu()`

## Purpose

Simulates a three-pulse DEER echo for a Cu(II)–NO two-electron system at X-band. The calculation performs brute-force time propagation and numerical powder averaging in Liouville space.

## Model and calculation

The two electron spins are separated by 20 Å and use distinct anisotropic Zeeman tensors. The script builds probe- and pump-spin excitation operators, evaluates the hard-echo sequence with delays `ta = 0.2 µs`, `tb = 0.1 µs`, and `tc = 25 ns`, and averages over a 1,600-point spherical grid. It plots the imaginary part of the echo against time in microseconds. The source estimates a run time of seconds.
