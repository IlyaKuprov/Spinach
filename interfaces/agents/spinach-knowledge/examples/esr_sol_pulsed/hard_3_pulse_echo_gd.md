# examples/esr_sol_pulsed/hard_3_pulse_echo_gd.m

- Signature: `hard_3_pulse_echo_gd()`

## Purpose

Calculates a three-pulse DEER echo for two Gd(III) centres. The source specifies hard pulses: the probe addresses the outermost zero-field-splitting (ZFS) transition and the pump the central transition. It warns that the Gd echo is very sharp and difficult to capture, and that experimental ZFS distributions are not included in this model.

## Spin system and fixed inputs

The field is 3.5 T. Both spins are E8 (eight-level electron spins) with scalar g value 2.002319. Each centre has an axial ZFS tensor with entries [D, D, −2D]/3, where D is entered as 0.56e9 in the example; the tensor Euler inputs are [0, 0, 0] and [0, π/2, 0]. The coordinates are [0, 0, 0] and 29.5 × [0, 1, 0]; the example does not state the coordinate unit. The full Zeeman-Hilbert basis is used (zeeman-hilb, approximation none). No parameter sweep is defined.

## Pulse sequence and sampled signal

The source embeds separate probe- and pump-pulse operators for the outermost and central E8 transitions, respectively. It prepares Lz magnetisation, detects with L+ on spin 1, and calls the hard-echo helper through powder averaging. The helper applies probe π/2, evolves for tb, pump π, evolves for ta−tb, probe π, and samples the centred echo window of width tc. The fixed delays are ta = 2 µs, tb = 1 µs, and tc = 50 ns, with 500 sampling steps; the spherical powder grid is rep_2ang_1600pts_sph.

The script plots the imaginary echo response against time in microseconds and creates no saved data file. It reports an estimated calculation time of seconds. Its stated sharp-echo and missing-ZFS-distribution caveat limits comparison with experimental lineshapes.

## Implementation

See the [example source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/hard_3_pulse_echo_gd.m) and the [three-pulse hard-echo helper](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_dipolar/deer_3p_hard_echo.m).
