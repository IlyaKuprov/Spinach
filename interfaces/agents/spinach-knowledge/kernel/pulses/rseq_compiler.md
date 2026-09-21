# kernel/pulses/rseq_compiler.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/rseq_compiler.m`
- Signature: `[P,T]=rseq_compiler(spin_system,L,Sx,Sy,pulse_phi,...`
- Total lines: 136

## Purpose

R sequence compiler. Uses the fact that R-sequences are very repetitive to pre-compile the minimal number of pulse propa- gators. Syntax: [P,T]=rseq_compiler(spin_system,L,Sx,Sy,pulse_phi,... pulse_amp,pulse_dur,element_type)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- L -background Liouvillian
- Sx,Sy -Cartesian spin operators pertaining to
- the spins affected by the pulses
- pulse_phi -the sequence of pulse phases, radians
- pulse_amp -RF nutation frequency in rad/s, a scalar
- because R-sequences are phase-modulated
- pulse_dur -duration of the pulses in the sequence
- element, a vector with the length mat-
- ching the number of pulses in the sequ-
- ence element (seconds)
- element_type -R element needs to be an inversion
- pulse; common ones are:
- '180_pulse' : simple inversion pulse
- '90270_pulse' : composite inversion pulse

## Outputs

- P -unique propagators, a cell array of matrices
- T -an index array of the same dimension as pulse_phi,
- specifying which propagator is to be used at which
- slice of the phase sequence

## Implementation structure

- R sequence compiler. Uses the fact that R-sequences are very
- repetitive to pre-compile the minimal number of pulse propa-
- gators. Syntax:
- [P,T]=rseq_compiler(spin_system,L,Sx,Sy,pulse_phi,...
- pulse_amp,pulse_dur,element_type)
- L -background Liouvillian
- Sx,Sy -Cartesian spin operators pertaining to
- the spins affected by the pulses
- pulse_phi -the sequence of pulse phases, radians
- pulse_amp -RF nutation frequency in rad/s, a scalar
- because R-sequences are phase-modulated
- pulse_dur -duration of the pulses in the sequence

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `phi()`, `propagator()`, `pulse_dur()`, `pulse_phi()`, `phi_dur()`, `ischar()`, `ishermitian()`, `isscalar()`, `any()`, `ismember()`.
