# kernel/pulses/rseq_compiler.m

- Signature: `[P,T]=rseq_compiler(spin_system,L,Sx,Sy,pulse_phi,pulse_amp,pulse_dur,element_type)`

## Purpose

Compiles the distinct pulse propagators needed by an R-sequence, reusing a propagator wherever the phase (and, for composite pulses, the phase-duration pair) repeats.

## Algorithm

For `'180_pulse'`, the function finds the unique values in `pulse_phi` and builds one propagator for each phase from `L+pulse_amp*(Sx*cos(phi)+Sy*sin(phi))`, using the single pulse duration. For `'90270_pulse'`, it pairs the alternating durations with the phases, finds unique phase-duration rows, and builds one propagator for each row.

## Parameters / inputs

- `spin_system` — Spinach system description used by `propagator`.
- `L` — background Liouvillian.
- `Sx`, `Sy` — Cartesian spin operators for the spins affected by the pulses; the source requires Hermitian matrices.
- `pulse_phi` — phase sequence, in radians.
- `pulse_amp` — scalar RF nutation frequency, in radians per second.
- `pulse_dur` — pulse duration in seconds; scalar for `'180_pulse'`, or two durations for `'90270_pulse'`.
- `element_type` — `'180_pulse'` for simple inversion pulses or `'90270_pulse'` for composite inversion pulses.

## Outputs

- `P` — cell array of unique propagator matrices.
- `T` — index array with the same dimensions as `pulse_phi`, indicating which entry of `P` to use at each phase-sequence slice.

[Spinach wiki page](https://spindynamics.org/wiki/index.php?title=rseq_compiler.m)
