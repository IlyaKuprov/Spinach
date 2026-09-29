# examples/esr_sol_pulsed/hard_3_pulse_echo_no.m

- Signature: `hard_3_pulse_echo_no()`

## Purpose

Calculates a three-pulse DEER spin echo for two nitroxide radicals at X-band, separated by 25 Å. The script uses brute-force time propagation and numerical powder averaging. Although its header describes averaging in Liouville space, the actual setup selects a Zeeman-Hilbert basis; the latter specifies the basis used by this example.

## Spin system and fixed inputs

The pair is modelled at 0.33 T with two electron spins (E) and identical anisotropic g-tensor principal values [2.0089, 2.0061, 2.0027]. The source assigns Euler values [1, 2, 3] and [3, 1, 2] to the tensors; it does not label their angle units. Coordinates are [0, 0, 0] and [25, 0, 0] Å. The exact Zeeman-Hilbert basis is selected with no basis approximation. These fixed source values define one pair, not a scan.

## Pulse sequence and sampled signal

The initial state is Lz magnetisation, and detection uses the probe-spin L+ coil state. Probe and pump operators are transverse Lx-type operators localised on spins 1 and 2. The script calls the hard-echo helper through powder averaging. Its ideal hard-pulse sequence is probe π/2, evolution for tb, pump π, evolution for ta−tb, probe π, then a centred echo window of width tc. The delays are ta = 1 µs, tb = 0.5 µs, and tc = 0.5 µs; it samples 256 steps and averages over rep_2ang_3200pts_sph.

The figure plots the imaginary echo response against a time axis from −tc/2 to +tc/2, in microseconds. It does not write a data file. The nitroxide g tensor is referenced to [DOI: 10.1063/1.1697233](http://dx.doi.org/10.1063/1.1697233); the example estimates a calculation time of seconds.

## Implementation

See the [example source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/hard_3_pulse_echo_no.m) and the [three-pulse hard-echo helper](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_dipolar/deer_3p_hard_echo.m).
