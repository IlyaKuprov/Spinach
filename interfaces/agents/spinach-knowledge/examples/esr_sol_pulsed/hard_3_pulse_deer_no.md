# examples/esr_sol_pulsed/hard_3_pulse_deer_no.m

- Signature: `hard_3_pulse_deer_no()`

## Purpose

Calculates a three-pulse DEER trace for two nitroxide radicals separated by 25 Å at X-band. The model uses brute-force time propagation and numerical powder averaging, including the orientation dependence of the dipolar interaction through the anisotropic g tensors.

## Model and calculation

The two electron spins are simulated at 0.33 T with the nitroxide g-tensor principal values and orientations specified in the source. The calculation uses a Zeeman Hilbert-space basis, a 3,200-point spherical powder grid, and 100 time steps of 10 ns. The resulting DEER trace is plotted against time in microseconds. The source estimates a run time of seconds.

Nitroxide g-tensor reference: http://dx.doi.org/10.1063/1.1697233
