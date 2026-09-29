# examples/esr_sol_pulsed/hard_3_pulse_echo_cu.m

- Signature: `hard_3_pulse_echo_cu()`

## Purpose

Calculates the three-pulse DEER spin echo for a Cu(II)–NO pair at X-band. The source describes brute-force time propagation with numerical powder averaging in Liouville space. This script calculates an echo window, not a scanned DEER distance distribution.

## Spin system and fixed inputs

The model contains two electron spins at 0.33 T. Their principal g values are [2.056, 2.056, 2.205] for the Cu(II) centre and [2.009, 2.006, 2.003] for the NO centre; both Euler inputs are [0, 0, 0]. Their coordinates are [0, 0, 0] and [20, 0, 0]; the example does not state the coordinate unit. It uses the full Zeeman-Liouville basis (sphten-liouv, approximation none); the source also disables trajectory-level SSR algorithms. The file fixes one field, geometry, tensor pair, and sequence rather than scanning them.

## Pulse sequence and sampled signal

The initial state is electron Lz; the coil detects L+ on spin 1. The probe and pump operators are transverse Lx-type operators localised on spins 1 and 2, respectively. The hard-echo helper, called under powder averaging, applies an ideal probe π/2, evolves for tb, applies a pump π, evolves for ta−tb, applies a probe π, then evolves to and samples a centred echo window of width tc. Here ta = 0.2 µs, tb = 0.1 µs, and tc = 25 ns; the sampled window uses 256 steps. The powder grid is rep_2ang_1600pts_sph.

The plotted observable is the imaginary echo response over −tc/2 to +tc/2, displayed in microseconds. The script creates a figure and does not specify a saved data file. Its comment estimates seconds for calculation time.

## Implementation

See the [example source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/esr_sol_pulsed/hard_3_pulse_echo_cu.m) and the [three-pulse hard-echo helper](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/esr_dipolar/deer_3p_hard_echo.m).
